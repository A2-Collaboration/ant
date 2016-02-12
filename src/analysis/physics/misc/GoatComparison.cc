#include "GoatComparison.h"

#include "root-addons/cbtaps_display/TH2CB.h"
#include "root-addons/cbtaps_display/TH2TAPS.h"


#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"

#include "TCanvas.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

using namespace std;

GoatComparison::GoatComparison(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    writeEvents(opts->Get<bool>("WriteEvents", false))
{
    steps = HistFac.makeTH1D("Steps","","",BinSettings(10),"steps");

    h_CBSumVetoE = HistFac.makeTH1D("CB Sum VetoE","E_{Veto} / MeV","",BinSettings(100,0,20),"h_CBSumVetoE");
    h_PIDSumE = HistFac.makeTH1D("PID Sum E","E / MeV","",BinSettings(100,0,20),"h_PIDSumE");


    n_photon_high = HistFac.make<TH2CB>("n_photon_high","n_photon_high");
    n_photon_low = HistFac.make<TH2CB>("n_photon_low","n_photon_low");

    IM_gg = HistFac.makeTH1D("IM(2#gamma)","IM [MeV]","",BinSettings(1000,0,1100),"IM_gg");

    BinSettings bins_photon_E(500,0,1000);


    photon_thetaE = HistFac.makeTH2D("Photon #theta vs. E","#theta / #circ","E / MeV",
                                     BinSettings(300,20,180),bins_photon_E,"photon_thetaE");

    BinSettings bins_clusterSize(20,0,20);

    photon_clusterSizeE = HistFac.makeTH2D("Photon ClusterSize vs. E","ClusterSize","E / MeV",
                                           bins_clusterSize,bins_photon_E,"photon_clusterSizeE");
}

GoatComparison::~GoatComparison()
{}

void GoatComparison::ProcessEvent(const TEvent& event, manager_t& manager)
{
    steps->Fill("Seen",1.0);

    const auto& cands = event.Reconstructed->Candidates;
    TCandidateList cands_cb;

    double CBSumVetoE = 0;
    for(const auto& c : cands) {
        if(c->Detector & Detector_t::Any_t::CB_Apparatus) {
            cands_cb.emplace_back(c);
            CBSumVetoE += c->VetoEnergy;
        }
    }

    double PIDSumE = 0;
    for(const TClusterPtr& cl : event.Reconstructed->Clusters) {
        if(cl->DetectorType == Detector_t::Type_t::PID)
            PIDSumE += cl->Energy;
    }

    const auto nCB = cands_cb.size();

    if(nCB != 2)
        return;

    steps->Fill("nCB==2",1.0);

    TCandidatePtr photon_high;
    TCandidatePtr photon_low;
    if(cands_cb.front()->CaloEnergy > cands_cb.back()->CaloEnergy) {
        photon_high = cands_cb.front();
        photon_low = cands_cb.back();
    }
    else {
        photon_high = cands_cb.back();
        photon_low = cands_cb.front();
    }

    TClusterPtr photon_high_cluster = photon_high->FindCaloCluster();
    TClusterPtr photon_low_cluster = photon_low->FindCaloCluster();

    if(!photon_high_cluster || !photon_low_cluster)
        return;
    steps->Fill("Clusters ok", 1.0);

    if(writeEvents)
        manager.SaveEvent();

    h_CBSumVetoE->Fill(CBSumVetoE);
    h_PIDSumE->Fill(PIDSumE);

    TLorentzVector sum(TParticle(ParticleTypeDatabase::Photon, photon_high)
                     + TParticle(ParticleTypeDatabase::Photon, photon_low));
    IM_gg->Fill(sum.M());

    n_photon_high->FillElement(photon_high_cluster->CentralElement, 1.0);
    n_photon_low->FillElement(photon_low_cluster->CentralElement, 1.0);

    photon_thetaE->Fill(std_ext::radian_to_degree(photon_high->Theta), photon_high->CaloEnergy);
    photon_thetaE->Fill(std_ext::radian_to_degree(photon_low->Theta), photon_low->CaloEnergy);

    photon_clusterSizeE->Fill(photon_high->ClusterSize, photon_high->CaloEnergy);
    photon_clusterSizeE->Fill(photon_low->ClusterSize, photon_low->CaloEnergy);

}

void GoatComparison::ShowResult()
{
    canvas()
            << IM_gg
            << drawoption("colz")
            << padoption::enable(padoption::LogZ)
            << n_photon_high
            << n_photon_low
            << photon_thetaE
            << photon_clusterSizeE
            << h_CBSumVetoE
            << h_PIDSumE
            << endc;
}

AUTO_REGISTER_PHYSICS(GoatComparison)
