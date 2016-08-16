#include "ReconstructCheck.h"

#include "base/Logger.h"

#include "TH1D.h"
#include "TH2D.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

ReconstructCheck::ReconstructCheck(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    Multiplicities(opts->Get<decltype(Multiplicities)>("Mult"))
{
    LOG(INFO) << "Multiplicities: " << Multiplicities;
    hMultiplicities = HistFac.makeTH1D("Multiplicity Distribution","#Candidates","%",BinSettings(10),"hMultiplicities");

    t.CreateBranches(HistFac.makeTTree("Tree"));
}




void ReconstructCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    t.Multiplicity = event.Reconstructed().Candidates.size();
    hMultiplicities->Fill(t.Multiplicity);

    t.Thetas().resize(t.Multiplicity);
    t.Energies().resize(t.Multiplicity);
    auto it_theta = t.Thetas().begin();
    auto it_energy = t.Energies().begin();

    for(const TCandidate& cand : event.Reconstructed().Candidates) {
        *it_theta = std_ext::radian_to_degree(cand.Theta);
        *it_energy = cand.CaloEnergy;
        ++it_theta;
        ++it_energy;
    }

    t.Tree->Fill();
}

void ReconstructCheck::Finish()
{
    hMultiplicities->Scale(100.0/hMultiplicities->GetEntries());
}

void ReconstructCheck::ShowResult()
{
    canvas(GetName()) << hMultiplicities << endc;
}


AUTO_REGISTER_PHYSICS(ReconstructCheck)
