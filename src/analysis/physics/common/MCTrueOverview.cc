#include "physics/common/MCTrueOverview.h"

#include "plot/root_draw.h"
#include "utils/particle_tools.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCTrueOverview::MCTrueOverview(const std::string& name, OptionsPtr opts):
    Physics(name,opts)
{
}

bool TryFindParticleTypeTree(const TParticleTree_t& ptree, ParticleTypeTreeDatabase::Channel& channel) {
    for(auto ch : ParticleTypeTreeDatabase()) {
        // channel is of ParticleTypeTreeDatabase::Channel
        auto ptree_db = ParticleTypeTreeDatabase::Get(ch);
        if(ptree->IsEqual(ptree_db, utils::ParticleTools::MatchByParticleName)) {
            channel = ch;
            return true;
        }
    }
    return false;
}

void MCTrueOverview::ProcessEvent(const TEvent& event, manager_t&)
{
    auto& ptree = event.MCTrue().ParticleTree;
    if(!ptree)
        return;
    ParticleTypeTreeDatabase::Channel channel;
    if(!TryFindParticleTypeTree(ptree, channel)) {
        LOG_N_TIMES(100, WARNING) << "Cannot find " << utils::ParticleTools::GetDecayString(ptree, false) << " in database (max 100x printed)";
        return;
    }


}

void MCTrueOverview::ShowResult()
{

}


AUTO_REGISTER_PHYSICS(MCTrueOverview)