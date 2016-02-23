#include "physics/common/CandidatesAnalysis.h"
#include "TLorentzVector.h"
#include <cmath>
#include "base/std_ext/math.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

CandidatesAnalysis::CandidatesAnalysis(const std::string& name, OptionsPtr opts):
    Physics(name,opts)
{
    nCandidatesEvent = HistFac.makeTH1D("Candidates/Event","Candidates","",BinSettings(15),"nCand");
    CandMultiplicities = HistFac.makeTH1D("Candidates Multi","","",BinSettings(15),"CandMult");

    for(int i=1; i<CandMultiplicities->GetNbinsX(); ++i) {
        CandMultiplicities->GetXaxis()->SetBinLabel(i,string(to_string(i-1)+"+").c_str());
    }

    energy = HistFac.makeTH1D("Energy","E [MeV]","",BinSettings(1000),"energy");
    theta  = HistFac.makeTH1D("Theta","#theta [#circ]","",BinSettings(360,0,180),"theta");
    phi    = HistFac.makeTH1D("Phi","#phi [#circ]","",BinSettings(720,-180,180),"phi");
    ggIM         = HistFac.makeTH1D("2 Neutral Candidates IM","M [MeV]","",BinSettings(1000),"ggIM");
    ttIM         = HistFac.makeTH1D("2 Candidates IM","M [MeV]","",BinSettings(1000),"ttIM");
    cbdEE = HistFac.makeTH2D("CB dE-E","E_{CB} [MeV]","dE_{PID} [MeV]", BinSettings(1000),BinSettings(100,0,30),"cb_dEE");
    cbtof = HistFac.makeTH2D("CB ToF","t_{CB} [ns]","E_{CB} [MeV]", BinSettings(300,-15,15),BinSettings(1000),"cb_tof");
    tapsdEE = HistFac.makeTH2D("TAPS dE-E","E_{TAPS} [MeV]","dE_{TAPSVeto} [MeV]", BinSettings(1000),BinSettings(100,0,30),"taps_dEE");
    tapstof = HistFac.makeTH2D("TAPS ToF","t_{TAPS} [ns]","E_{TAPS} [MeV]", BinSettings(300,-15,15),BinSettings(1000),"taps_tof");
    detectors = HistFac.makeTH1D("Detectors","","", BinSettings(1),"detectors");
    psa = HistFac.makeTH2D("TAPS PSA (Charged)","E_{long} [MeV]","E_{short} [MeV]", BinSettings(1000),BinSettings(1000),"psa");
    psa_all = HistFac.makeTH2D("TAPS PSA","E_{long} [MeV]","E_{short} [MeV]", BinSettings(1000),BinSettings(1000),"psa_all");
    psa_all_angles = HistFac.makeTH2D("TAPS PSA","#phi [#circ]","r", BinSettings(160,45-20,45+20),BinSettings(250,0,500),"psa_all_angles");
}

void CandidatesAnalysis::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& candidates = event.Reconstructed->Candidates;

    nCandidatesEvent->Fill(candidates.size());

    for(size_t i=0; i<=candidates.size(); ++i) {
        CandMultiplicities->Fill(i);
    }

    auto i = candidates.begin();

    while(i!=candidates.end()) {
        const TCandidatePtr& ci = *i;
        energy->Fill((*i)->CaloEnergy);
        theta->Fill((*i)->Theta*TMath::RadToDeg());
        phi->Fill((*i)->Phi*TMath::RadToDeg());
        detectors->Fill(string(ci->Detector).c_str(),1);

        if((*i)->CaloEnergy>20.0) {

            if(ci->Detector & Detector_t::Any_t::CB_Apparatus) {
                cbdEE->Fill(ci->CaloEnergy,ci->VetoEnergy);
                cbtof->Fill(ci->Time, ci->CaloEnergy);
            } else if(ci->Detector & Detector_t::Any_t::TAPS_Apparatus) {
                tapsdEE->Fill(ci->CaloEnergy,ci->VetoEnergy);
                tapstof->Fill(ci->Time, ci->CaloEnergy);


                // extract short gate stuff

                const auto& cluster = ci->FindCaloCluster();

                if(cluster)

                    psa_all_angles->Fill(std_ext::radian_to_degree(cluster->GetPSAAngle()), cluster->GetPSARadius());

                    for(const TClusterHit& clusterhit : cluster->Hits) {

                        if(clusterhit.Channel == cluster->CentralElement) {
                            double central_e = 0.0;
                            for(const TClusterHit::Datum& datum : clusterhit.Data) {

                                if(datum.Type == Channel_t::Type_t::Integral)
                                    central_e = datum.Value;

                                if(datum.Type == Channel_t::Type_t::IntegralShort) {

                                    if(ci->VetoEnergy<0.5)
                                        psa->Fill(central_e, datum.Value);

                                    psa_all->Fill(central_e, datum.Value);

                                }
                            }
                        }
                    }
            }

            const TParticle a(ParticleTypeDatabase::Photon,*i);
            auto j = i;
            ++j;
            while(j!=candidates.end()) {
                if((*i)->CaloEnergy>20.0) {
                    const TParticle b(ParticleTypeDatabase::Photon,*j);
                    const TLorentzVector s = a + b;

                    if((*i)->VetoEnergy==0 && (*j)->VetoEnergy==0)
                        ggIM->Fill(s.M());

                    ttIM->Fill(s.M());
                }
                ++j;
            }
        }
        ++i;
    }
}


void CandidatesAnalysis::Finish()
{

}


void CandidatesAnalysis::ShowResult()
{
    canvas("CandidatesAnalysis")
            << ggIM << energy << theta << phi
            << nCandidatesEvent << CandMultiplicities
            << detectors
            << drawoption("colz")
            << cbdEE << cbtof
            << tapsdEE << tapstof
            << psa << psa_all << psa_all_angles
            << endc;
}

AUTO_REGISTER_PHYSICS(CandidatesAnalysis)
