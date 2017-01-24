#include "TestMCWeighting.h"

#include "utils/MCWeighting.h"
#include "base/Logger.h"

#include "TDirectory.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

namespace ant {
namespace analysis {
namespace utils {
struct MCWeighting_test :  MCWeighting {
    // use ctor from base class
    using MCWeighting::MCWeighting;
    // expose some actually protected methods
    // for testing/verification here...
    using MCWeighting::GetBeamE;
    using MCWeighting::GetCosTheta;

};
}}} // namespace ant::analysis::utils

TestMCWeighting::TestMCWeighting(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    mcWeighting(std_ext::make_unique<utils::MCWeighting_test>(
                    HistFac,
                    utils::MCWeighting::EtaPrime)
                )
{
    t.CreateBranches(HistFac.makeTTree("t"));
}

void TestMCWeighting::ProcessEvent(const TEvent& event, manager_t&)
{
    auto& tree = event.MCTrue().ParticleTree;
    t.BeamE = mcWeighting->GetBeamE(tree);
    t.CosTheta = mcWeighting->GetCosTheta(tree);
    t.Tree->Fill();

    mcWeighting->SetParticleTree(tree);
    mcWeighting->Fill();
}

void TestMCWeighting::ShowResult()
{
    mcWeighting->FriendTTree(t.Tree);

    ant::canvas(GetName())
            << TTree_drawable(t.Tree, "BeamE")
            << TTree_drawable(t.Tree, "CosTheta")
            << drawoption("colz")
            << TTree_drawable(t.Tree, "CosTheta:BeamE >> h","MCWeight")
            << endc;

    // ugly way of getting our histogram from TTree_drawable
    auto h = dynamic_cast<TH2F*>(gDirectory->Get("h"));
    LOG(INFO) << "Total integral of 2D plotted weigths (should equal number of input events): "
              << h->Integral(1,h->GetXaxis()->GetNbins()-1,
                             1,h->GetYaxis()->GetNbins()-1);
}

void TestMCWeighting::Finish()
{
    mcWeighting->Finish();
}

AUTO_REGISTER_PHYSICS(TestMCWeighting);