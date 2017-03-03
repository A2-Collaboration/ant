#include "physics/common/MCClusteringCheck.h"

#include "plot/root_draw.h"
#include "utils/particle_tools.h"
#include "base/Logger.h"

#include "TH1D.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCClusteringCheck::MCClusteringCheck(const std::string& name, OptionsPtr opts):
    Physics(name,opts),
    opening_angles(
        [this] () {
    std::remove_const<decltype(opening_angles)>::type v;
    v.push_back({{10, 15}, HistFac});
    v.push_back({{15, 20}, HistFac});
    v.push_back({{20, 25}, HistFac});
    v.push_back({{25, 30}, HistFac});
    v.push_back({{30, 50}, HistFac});
    return v;
}())
{
    h_Steps = HistFac.makeTH1D("Steps","","",BinSettings(10),"h_Steps");
}


void MCClusteringCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& true_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
    const auto& true_photons = true_particles.Get(ParticleTypeDatabase::Photon);
    if(true_photons.size()<2)
        return;

    h_Steps->Fill("TruePhotons>=2",1.0);
    h_Steps->Fill("TruePhotons==2",true_photons.size()==2);

    auto true_photon1 = true_photons.front();
    auto true_photon2 = true_photons.back();
    if(true_photon1->Ek() < true_photon2->Ek())
        std::swap(true_photon1, true_photon2);

    const auto opening_angle = std_ext::radian_to_degree(true_photon1->Angle(*true_photon2));

    const auto& cands = event.Reconstructed().Candidates;
    const auto nCands = cands.size();

    // hm, maybe not the best matching procedure
    struct matched_candidate_t {
        explicit matched_candidate_t(const TParticle& p) : true_photon(p) {}
        const TParticle& true_photon;
        TCandidatePtr cand;
        double min_angle = std_ext::inf;
        void test(const TCandidatePtr& c) {
            auto angle = true_photon.Angle(*c);
            // take it if it has smaller opening angle
            if(angle < min_angle) {
                min_angle = angle;
                cand = c;
            }
        }
    };

    matched_candidate_t best_cand1(*true_photon1);
    matched_candidate_t best_cand2(*true_photon2);

    for(auto cand : cands.get_iter()) {
        if(cand->Detector & Detector_t::Type_t::CB) {
            best_cand1.test(cand);
            best_cand2.test(cand);
        }
    }

    if(best_cand1.cand == best_cand2.cand) {
        h_Steps->Fill("Cand1==Cand2",1.0);
        best_cand1.cand = nullptr;
        best_cand2.cand = nullptr;
    }

    const double maximum_match_angle = std_ext::degree_to_radian(15.0);
    if(best_cand1.min_angle > maximum_match_angle) {
        h_Steps->Fill("Cand1 too far",1.0);
        best_cand1.cand = nullptr;
    }
    if(best_cand2.min_angle > maximum_match_angle) {
        h_Steps->Fill("Cand2 too far",1.0);
        best_cand2.cand = nullptr;
    }

    if(best_cand1.cand && best_cand2.cand) {
        h_Steps->Fill("Cands match", 1.0);
    }

    for(auto& item : opening_angles) {
        // stop filling once fitting opening angle is found
        if(item.Fill(opening_angle, nCands, *true_photon1, *true_photon2,
                     best_cand1.cand, best_cand2.cand))
            break;
    }
}


void MCClusteringCheck::Finish()
{

}

void MCClusteringCheck::ShowResult()
{
    canvas c(GetName()+": Opening Angle Bins");
    for(auto& item : opening_angles) {
        item.Show(c);
        c << endr;
    }
    c << endc;
    canvas(GetName()+": Overview") << h_Steps << endc;
}


MCClusteringCheck::opening_angle_t::opening_angle_t(const interval<double> opening_angle_range_,
                                                    const HistogramFactory& HistFac) :
    opening_angle_range(opening_angle_range_)
{
    string name = std_ext::formatter() << "OpAng" << opening_angle_range;
    HistogramFactory histFac(name, HistFac, name); // use name as titleprefix

    const BinSettings bins_Theta(30, 20, 160);
    const BinSettings bins_EtrueErec(50, 0.5, 1.3);
    const BinSettings bins_OpeningAngle(50, 0, 20);

    h_nCands = histFac.makeTH1D("nCands", "nCands", "", BinSettings(5), "nCands");

    h_ErecEtrue1 = histFac.makeTH2D("E_{rec}/E_{true} 1","#theta_{true} / #circ","E_{rec}/E_{true}",
                                    bins_Theta, bins_EtrueErec, "h_ErecEtrue1");
    h_ErecEtrue2 = histFac.makeTH2D("E_{rec}/E_{true} 2","#theta_{true} / #circ","E_{rec}/E_{true}",
                                    bins_Theta, bins_EtrueErec, "h_ErecEtrue2");

    h_OpeningAngle1 = histFac.makeTH2D("OpeningAngle 1","#theta_{true} / #circ","Opening Angle / #circ",
                                       bins_Theta, bins_OpeningAngle, "h_OpeningAngle1");
    h_OpeningAngle2 = histFac.makeTH2D("OpeningAngle 2","#theta_{true} / #circ","Opening Angle / #circ",
                                       bins_Theta, bins_OpeningAngle, "h_OpeningAngle2");
}

bool MCClusteringCheck::opening_angle_t::Fill(
        double opening_angle, unsigned nCands,
        const TParticle& true_photon1, const TParticle& true_photon2,
        const TCandidatePtr& best_cand1, const TCandidatePtr& best_cand2) const
{
    if(!opening_angle_range.Contains(opening_angle))
        return false;

    h_nCands->Fill(nCands);

    if(best_cand1 && best_cand2) {
        const auto true_Theta1 = std_ext::radian_to_degree(true_photon1.Theta());
        const auto true_Theta2 = std_ext::radian_to_degree(true_photon2.Theta());
        const auto true_Ek1 = true_photon1.Ek();
        const auto true_Ek2 = true_photon2.Ek();

        const auto rec_Ek1 = best_cand1->CaloEnergy;
        const auto rec_Ek2 = best_cand2->CaloEnergy;

        h_ErecEtrue1->Fill(true_Theta1, rec_Ek1/true_Ek1);
        h_ErecEtrue2->Fill(true_Theta2, rec_Ek2/true_Ek2);

        h_OpeningAngle1->Fill(true_Theta1, std_ext::radian_to_degree(true_photon1.Angle(*best_cand1)));
        h_OpeningAngle2->Fill(true_Theta2, std_ext::radian_to_degree(true_photon2.Angle(*best_cand2)));

    }

    return true;
}

void MCClusteringCheck::opening_angle_t::Show(canvas& c) const
{
    c << h_nCands
      << drawoption("colz")
      << h_ErecEtrue1 << h_ErecEtrue2
      << h_OpeningAngle1 << h_OpeningAngle2;
}



AUTO_REGISTER_PHYSICS(MCClusteringCheck)
