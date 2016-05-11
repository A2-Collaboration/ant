#include "pi0eta_lost_g.h"
#include "base/Logger.h"
#include <algorithm>
#include "base/std_ext/math.h"

#include "TTree.h"

#include "utils/matcher.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;

Pi0EtaLostG::Pi0EtaLostG(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    tree(HistFac.makeTTree("tree"))
{
    t.CreateBranches(tree);
}

Pi0EtaLostG::~Pi0EtaLostG()
{
}

void Pi0EtaLostG::ProcessEvent(const TEvent& event, manager_t& manager)
{

    const auto& mc_cands = event.MCTrue().Particles.GetAll();
    const auto& re_cands = event.Reconstructed().Candidates;

    if(mc_cands.size() != 5 || re_cands.size() != 4)
        return;

    // find matched and unmatched
}



Pi0EtaLostG::Tree_t::Tree_t()
{

}

AUTO_REGISTER_PHYSICS(Pi0EtaLostG)
