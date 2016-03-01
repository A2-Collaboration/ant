#include "PID_Energy_etaDalitz.h"

#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

template<typename T>
void PID_Energy_etaDalitz::shift_right(std::vector<T>& v)
{
    std::rotate(v.begin(), v.end() -1, v.end());
}

PID_Energy_etaDalitz::PID_Energy_etaDalitz(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);

    const BinSettings cb_channels(detector->GetNChannels());
    const BinSettings energybins(1000, 0, 10);

    eegPID = HistFac.makeTH2D("2 charged 1 neutral (PID,PID)", "PID Energy [MeV]", "#",
                              energybins, cb_channels, "eegPID");
    etaIM = HistFac.makeTH1D("#eta IM all comb", "IM [MeV]", "#", BinSettings(1200), "etaIM");
    hCopl = HistFac.makeTH1D("coplanarity #eta - proton all comb", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), "hCopl");
    protonVeto = HistFac.makeTH1D("Veto energy identified proton", "Veto [MeV]", "#", energybins, "protonVeto");
}

void PID_Energy_etaDalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed->Candidates;

    if(cands.size() != 4)
        return;

    TLorentzVector eta;
    TLorentzVector proton;
    const interval<double> eta_im({ETA_IM-ETA_SIGMA, ETA_IM+ETA_SIGMA});
    const interval<double> coplanarity({-20, 20});
    vector<TCandidatePtr> comb;
    for( auto p : cands )
        comb.emplace_back(p);
    bool found;
    for( size_t j = 0; j < cands.size(); j++ ) {  // loop to test all different combinations
        found = true;
        eta.SetXYZT(0,0,0,0);
        for( size_t i = 0; i < comb.size()-1; i++ )
            eta += TParticle(ParticleTypeDatabase::Photon, comb.at(i));
        proton = TParticle(ParticleTypeDatabase::Proton, comb.back());
        const double copl = std_ext::radian_to_degree(abs(eta.Phi() - proton.Phi())) - 180.;
        etaIM->Fill(eta.M());
        hCopl->Fill(copl);
        if (eta_im.Contains(eta.M()) && coplanarity.Contains(copl)
                && std::count_if(comb.begin(), comb.end()-1, [](TCandidatePtr c){ return c->VetoEnergy; }) >= 2)
            break;
        found = false;
        shift_right(comb);
    }
    if(!found)
        return;

    protonVeto->Fill(comb.back()->VetoEnergy);
    // at this point a possible eta Dalitz candidate was found, work only with eta final state
    comb.pop_back();

    sort(comb.begin(), comb.end(),
         [] (const TCandidatePtr& a, const TCandidatePtr& b) {
            return a->VetoEnergy > b->VetoEnergy;
         });

    const TCandidatePtr& l1 = comb.at(0);
    const TCandidatePtr& l2 = comb.at(1);
    // suppress conversion decays
    if(l1->FindVetoCluster()->CentralElement == l2->FindVetoCluster()->CentralElement)
        return;
    const double eeIM = (TParticle(ParticleTypeDatabase::eMinus, l1) + TParticle(ParticleTypeDatabase::eMinus, l2)).M();
    // suppress pi0
    if(eeIM > 135.)
        return;

    for( const TCandidatePtr& c : comb )
        if(c->VetoEnergy)
            eegPID->Fill(c->VetoEnergy, c->FindVetoCluster()->CentralElement);
}

void PID_Energy_etaDalitz::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << eegPID << endc;
}

AUTO_REGISTER_PHYSICS(PID_Energy_etaDalitz)
