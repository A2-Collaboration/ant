#include "TestMCWeighting.h"

#include "utils/MCWeighting.h"

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
    mcWeighting(std_ext::make_unique<utils::MCWeighting_test>(HistFac, utils::MCWeighting::EtaPrime))
{
    t.CreateBranches(HistFac.makeTTree("t"));
}

void TestMCWeighting::ProcessEvent(const TEvent& event, manager_t&)
{
    auto& tree = event.MCTrue().ParticleTree;
    t.BeamE = mcWeighting->GetBeamE(tree);
    t.CosTheta = mcWeighting->GetCosTheta(tree);
    t.Tree->Fill();
}

void TestMCWeighting::ShowResult()
{
    ant::canvas(GetName())
            << TTree_drawable(t.Tree, "BeamE")
            << TTree_drawable(t.Tree, "CosTheta")
            << endc;
}

AUTO_REGISTER_PHYSICS(TestMCWeighting);