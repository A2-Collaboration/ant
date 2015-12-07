#include "ExtractResolutions.h"

#include "TTree.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

physics::ExtractResolutions::ExtractResolutions(const std::string& name, PhysOptPtr opts):
    Physics(name, opts)
{
    unsigned bins = 0;
    if(opts->Get<string>("Detector") == "CB") {
        det = Detector_t::Type_t::CB;
        bins = 720;
        HistFac.SetTitlePrefix("CB");
    } else if(opts->Get<string>("Detector") == "TAPS") {
        det = Detector_t::Type_t::TAPS;
        bins = 438;
        HistFac.SetTitlePrefix("TAPS");
    }
    const BinSettings b(bins);
    const BinSettings eb(1600);

    E_offset = HistFac.makeTH2D("E offset","Element", "E [MeV]", b, eb, "E_offset");
    E_sigma  = HistFac.makeTH2D("E #sigma","Element", "E [MeV]", b, eb, "E_sigma");

    Theta_offset = HistFac.makeTH2D("Theta offset","Element","E [MeV]", b, eb, "Theta_offset");
    Theta_sigma  = HistFac.makeTH2D("Theta sigma", "Element","E [MeV]", b, eb, "Theta_sigma");

    Phi_offset = HistFac.makeTH2D("Phi offset","Element","E [MeV]", b, eb, "Phi_offset");
    Phi_sigma  = HistFac.makeTH2D("Phi sigma", "Element","E [MeV]", b, eb, "Phi_sigma");

    tree = HistFac.makeTTree("tree");
    tree->Branch("E",     &b_E);
    tree->Branch("Theta", &b_Theta);
    tree->Branch("Phi",   &b_Phi);
    tree->Branch("Element", &b_Element);
}

void physics::ExtractResolutions::ProcessEvent(const data::Event& event)
{
    const auto& mcparticles = event.MCTrue.Particles.GetAll();

    if(mcparticles.size() == 1) {

        const auto& recparticles = event.Reconstructed.Particles.GetAll();

        if(recparticles.size() ==1 ) {

            const auto& mcp = mcparticles.front();
            const auto& rep = recparticles.front();

            if(rep->Candidate && rep->Candidate->Detector == det) {
                const auto& c = rep->Candidate->FindCaloCluster();
                if(c) {
                    b_E     = rep->Ek() - mcp->Ek();
                    b_Theta = rep->Theta() - mcp->Theta();
                    b_Phi   = rep->Phi() - mcp->Phi();
                    b_Element = c->CentralElement;
                    tree->Fill();
                }
            }
        }
    }
}

void physics::ExtractResolutions::ShowResult()
{

}

AUTO_REGISTER_PHYSICS(ExtractResolutions)
