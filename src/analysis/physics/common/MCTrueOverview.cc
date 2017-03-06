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

bool TryFindParticleTypeTree(const TParticleTree_t& ptree,
                             ParticleTypeTreeDatabase::Channel& channel,
                             ParticleTypeTree& typetree) {
    for(auto ch : ParticleTypeTreeDatabase()) {
        // channel is of ParticleTypeTreeDatabase::Channel
        auto db_typetree = ParticleTypeTreeDatabase::Get(ch);
        if(ptree->IsEqual(db_typetree, utils::ParticleTools::MatchByParticleName)) {
            channel = ch;
            typetree = db_typetree;
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
    ParticleTypeTree typetree;
    if(!TryFindParticleTypeTree(ptree, channel, typetree)) {
        LOG_N_TIMES(100, WARNING) << "Cannot find " << utils::ParticleTools::GetDecayString(ptree, false) << " in database (max 100x printed)";
        return;
    }

    // search for channel specific histograms
    auto it_perChannel = channels.find(channel);
    if(it_perChannel == channels.end()) {
        auto it = channels.emplace(make_pair(channel, perChannel_t(HistFac, typetree)));
        it_perChannel = it.first;
    }
    auto& perChannel = it_perChannel->second;
    perChannel.Fill(ptree, typetree);

}

void MCTrueOverview::ShowResult()
{

}

MCTrueOverview::perChannel_t::perChannel_t(const HistogramFactory& histFac, const ParticleTypeTree& typetree)
{
    const auto& decaystring = utils::ParticleTools::GetDecayString(typetree, true);
    HistogramFactory HistFac(decaystring, histFac, decaystring);

    histtree = typetree->DeepCopy<histnode_t>([HistFac] (const ParticleTypeTree& t) {
        // only create histograms at this node if a leaf daughter exists
        // indicated by non-empty leafTypes
        std::vector<histnode_t::typeptr_t> leafTypes;
        for(auto d  : t->Daughters()) {
            if(d->IsLeaf()) {
                leafTypes.push_back(addressof(d->Get()));
            }
        }
        auto subdecaystring = t->Get().PrintName();
        auto histFacPtr = leafTypes.empty() ? nullptr : std_ext::make_unique<const HistogramFactory>(subdecaystring, HistFac, subdecaystring);
        return histnode_t(move(histFacPtr), leafTypes);
    });

}

void MCTrueOverview::perChannel_t::Fill(const TParticleTree_t& ptree, const ParticleTypeTree& typetree) const
{

}

MCTrueOverview::perChannel_t::histnode_t::histnode_t(std::unique_ptr<const HistogramFactory> histFacPtr,
                                                     const vector<typeptr_t>& leafTypes)
{
    if(!histFacPtr)
        return;
    auto& HistFac = *histFacPtr;
    for(auto typeptr : leafTypes) {
        auto typestr = typeptr->PrintName();
        hists.emplace(make_pair(typeptr, HistogramFactory(typestr, HistFac, typestr)));
    }
}



MCTrueOverview::perChannel_t::histnode_t::perType_t::perType_t(const HistogramFactory& HistFac)
{
    const AxisSettings axis_Theta("#theta / #circ", {50, 0, 180});
    const AxisSettings axis_Ek("E_{k} / MeV", {100, 0, 1000});

    h_EkTheta = HistFac.makeTH2D("E_{k} vs. #theta", axis_Theta, axis_Ek, "h_EkTheta");
}

AUTO_REGISTER_PHYSICS(MCTrueOverview)