#include "MCClusteringCheck.h"

#include "utils/ParticleTools.h"
#include "plot/HistStyle.h"
#include "base/Logger.h"

#include "TH1D.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCClusteringCheck::MCClusteringCheck(const std::string& name, OptionsPtr opts):
    Physics(name,opts),
    detectorType(opts->Get<bool>("UseTAPS", false) ?
                     Detector_t::Type_t::TAPS : Detector_t::Type_t::CB),
    opening_angles(
        [this] () {
    std::remove_const<decltype(opening_angles)>::type v;
    v.push_back({{10, 15}, HistFac, detectorType});
    v.push_back({{15, 20}, HistFac, detectorType});
    v.push_back({{20, 25}, HistFac, detectorType});
    v.push_back({{25, 30}, HistFac, detectorType});
    v.push_back({{30, 50}, HistFac, detectorType});
    return v;
}())
{
    h_Steps = HistFac.makeTH1D("Steps",{"",BinSettings(10)},"h_Steps");
    h_Cands_OpAng.resize(5);
    for(int i=0;i<int(h_Cands_OpAng.size());++i) {
        h_Cands_OpAng[i] = HistFac.makeTH1D("nCand="+to_string(i-1),"#Delta#alpha / #circ","Fraction of total",BinSettings(49,1,50),"h_Cands_OpAng_"+to_string(i));
        h_Cands_OpAng[i]->SetLineColor(plot::histstyle::color_t::GetDark(i));
        h_Cands_OpAng[i]->SetLineWidth(2);
    }
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

    h_Cands_OpAng[0]->Fill(opening_angle);
    if(nCands<h_Cands_OpAng.size()-1)
        h_Cands_OpAng[nCands+1]->Fill(opening_angle);

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
        if(cand->Detector & detectorType) {
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

    if(opening_angle>10) {
        h_Steps->Fill("OpAng > 10#circ", 1.0);
    }
    if(best_cand1.cand && best_cand2.cand) {
        h_Steps->Fill("Cands match", 1.0);
        if(opening_angle>10) {
            h_Steps->Fill("Cands match > 10#circ", 1.0);
        }
    }

    TCandidatePtrList unmatched_cands;
    for(auto cand : cands.get_iter()) {
        if(cand->Detector & detectorType) {
            // cand is not assigned to anything
            if(best_cand1.cand != cand &&
               best_cand2.cand != cand)
                // then add it to unmatched
                unmatched_cands.emplace_back(cand);
        }
    }

    if(!unmatched_cands.empty())
        h_Steps->Fill("Unmatched cands>0",1.0);

    for(auto& item : opening_angles) {
        // stop filling once fitting opening angle is found
        if(item.Fill(opening_angle, nCands, *true_photon1, *true_photon2,
                     best_cand1.cand, best_cand2.cand, unmatched_cands))
            break;
    }
}

void MCClusteringCheck::ShowResult()
{
    canvas c(GetName()+": Opening Angle Bins");
    for(auto& item : opening_angles) {
        item.Show(c);
        c << endr;
    }
    c << endc;
    canvas c_overview(GetName()+": Overview");
    c_overview << h_Steps;

    for(int i=1;i<int(h_Cands_OpAng.size());++i) {
        c_overview << h_Cands_OpAng[i];
    }

    c_overview << endc;
}

void MCClusteringCheck::Finish()
{
    for(int i=1;i<int(h_Cands_OpAng.size());++i) {
         h_Cands_OpAng[i]->Divide(h_Cands_OpAng[0]);
    }
}

MCClusteringCheck::opening_angle_t::opening_angle_t(const interval<double> opening_angle_range_,
                                                    const HistogramFactory& HistFac,
                                                    Detector_t::Type_t detectorType) :
    opening_angle_range(opening_angle_range_)
{
    string name = std_ext::formatter() << "OpAng" << opening_angle_range;
    HistogramFactory histFac(name, HistFac, name); // use name as titleprefix

    const AxisSettings axis_TrueTheta("#theta_{true} / #circ",
                                      detectorType == Detector_t::Type_t::CB ?
                                          BinSettings{30, 20, 160} : BinSettings{30, 0, 25});
    const AxisSettings axis_TrueEnergy("E^{kin}_{true} / MeV", BinSettings{160, 0, 1600});
    const AxisSettings axis_EtrueErec("E_{rec}/E_{true}", {50, 0.7, 1.3});
    const AxisSettings axis_OpeningAngle("Opening Angle / #circ",
                                         detectorType == Detector_t::Type_t::CB ?
                                             BinSettings{50, 0, 12} : BinSettings{50, 0, 3});
    const AxisSettings axis_DiffAngleTheta("#theta_{rec} - #theta_{true} / #circ",
                                           detectorType == Detector_t::Type_t::CB ?
                                               BinSettings{50, -6, 6} : BinSettings{50, -3, 3});
    const AxisSettings axis_DiffAnglePhi("#phi_{rec} - #phi_{true} / #circ",
                                         detectorType == Detector_t::Type_t::CB ?
                                             BinSettings{50, -6, 6} : BinSettings{50, -3, 3});

    h_nCands = histFac.makeTH1D("nCands", {"nCands", BinSettings(8)}, "nCands");
    h_nSplits = histFac.makeTH1D("nSplits", {"nSplits", BinSettings(3)}, "nSplits");

    h_ErecEtrue1 = histFac.makeTH2D("E_{rec}/E_{true} 1", axis_TrueTheta, axis_EtrueErec, "h_ErecEtrue1");
    h_ErecEtrue2 = histFac.makeTH2D("E_{rec}/E_{true} 2", axis_TrueTheta, axis_EtrueErec, "h_ErecEtrue2");

    h_OpeningAngle1 = histFac.makeTH2D("OpAngle 1", axis_TrueTheta, axis_OpeningAngle, "h_OpeningAngle1");
    h_OpeningAngle2 = histFac.makeTH2D("OpAngle 2", axis_TrueTheta, axis_OpeningAngle, "h_OpeningAngle2");

    h_DiffAngleTheta1 = histFac.makeTH2D("DiffAngleTheta 1", axis_TrueEnergy, axis_DiffAngleTheta, "h_DiffAngleTheta1");
    h_DiffAngleTheta2 = histFac.makeTH2D("DiffAngleTheta 2", axis_TrueEnergy, axis_DiffAngleTheta, "h_DiffAngleTheta2");
    h_DiffAnglePhi1 = histFac.makeTH2D("DiffAnglePhi 1", axis_TrueEnergy, axis_DiffAnglePhi, "h_DiffAnglePhi1");
    h_DiffAnglePhi2 = histFac.makeTH2D("DiffAnglePhi 2", axis_TrueEnergy, axis_DiffAnglePhi, "h_DiffAnglePhi2");

    h_nUnmatchedCandsMinAngle = histFac.makeTH1D("nUnmatchedCandsMinAngle",
                                                 {"Min Angle / #circ",  detectorType == Detector_t::Type_t::CB ?
                                                  BinSettings{50,0,180} : BinSettings{50,0,50}},
                                                 "h_nUnmatchedCandsMinAngle");
}

bool MCClusteringCheck::opening_angle_t::Fill(
        double opening_angle, unsigned nCands,
        const TParticle& true_photon1, const TParticle& true_photon2,
        const TCandidatePtr& best_cand1, const TCandidatePtr& best_cand2,
        const TCandidatePtrList& unmatched_cands) const
{
    if(!opening_angle_range.Contains(opening_angle))
        return false;

    h_nCands->Fill(nCands);

    if(best_cand1 && best_cand2) {
        const auto true_Theta1 = std_ext::radian_to_degree(true_photon1.Theta());
        const auto true_Theta2 = std_ext::radian_to_degree(true_photon2.Theta());
        const auto true_Phi1 = std_ext::radian_to_degree(true_photon1.Phi());
        const auto true_Phi2 = std_ext::radian_to_degree(true_photon2.Phi());
        const auto true_Ek1 = true_photon1.Ek();
        const auto true_Ek2 = true_photon2.Ek();

        const auto rec_Ek1 = best_cand1->CaloEnergy;
        const auto rec_Ek2 = best_cand2->CaloEnergy;

        h_nSplits->Fill(best_cand1->FindCaloCluster()->HasFlag(TCluster::Flags_t::Split) +
                        best_cand2->FindCaloCluster()->HasFlag(TCluster::Flags_t::Split));

        h_ErecEtrue1->Fill(true_Theta1, rec_Ek1/true_Ek1);
        h_ErecEtrue2->Fill(true_Theta2, rec_Ek2/true_Ek2);

        h_OpeningAngle1->Fill(true_Theta1, std_ext::radian_to_degree(true_photon1.Angle(*best_cand1)));
        h_OpeningAngle2->Fill(true_Theta2, std_ext::radian_to_degree(true_photon2.Angle(*best_cand2)));

        h_DiffAngleTheta1->Fill(true_Ek1, std_ext::radian_to_degree(best_cand1->Theta) - true_Theta1);
        h_DiffAngleTheta2->Fill(true_Ek2, std_ext::radian_to_degree(best_cand2->Theta) - true_Theta2);

        h_DiffAnglePhi1->Fill(true_Ek1, std_ext::radian_to_degree(best_cand1->Phi) - true_Phi1);
        h_DiffAnglePhi2->Fill(true_Ek2, std_ext::radian_to_degree(best_cand2->Phi) - true_Phi2);

        // for each unmatched cand, fill the minimum angle to one of the true photons
        for(auto& cand : unmatched_cands) {
            const auto angle1 = std_ext::radian_to_degree(true_photon1.Angle(*cand));
            const auto angle2 = std_ext::radian_to_degree(true_photon2.Angle(*cand));
            h_nUnmatchedCandsMinAngle->Fill(angle1 < angle2 ? angle1 : angle2);
        }

    }

    return true;
}

void MCClusteringCheck::opening_angle_t::Show(canvas& c) const
{
    c << padoption::LogY << h_nCands
      << padoption::LogY  << h_nSplits
      << drawoption("colz")
      << h_ErecEtrue1 << h_ErecEtrue2
      << h_OpeningAngle1 << h_OpeningAngle2
      << h_DiffAngleTheta1 << h_DiffAngleTheta2
      << h_DiffAnglePhi1 << h_DiffAnglePhi2
      << padoption::LogY << h_nUnmatchedCandsMinAngle;
}



AUTO_REGISTER_PHYSICS(MCClusteringCheck)
