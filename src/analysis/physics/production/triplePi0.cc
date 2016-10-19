#include "triplePi0.h"

//#include "base/ParticleType.h"
//#include "base/ParticleTypeTree.h"
//#include "base/std_ext/math.h"

#include "expconfig/ExpConfig.h"

//#include "plot/root_draw.h"

//#include "utils/combinatorics.h"
#include "utils/particle_tools.h"
//#include "utils/ParticleID.h"

//#include <algorithm>
//#include <cassert>
//#include <chrono>


using namespace std;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const triplePi0::named_channel_t triplePi0::signal =
    {"3Pi0Prod",    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)};
const triplePi0::named_channel_t triplePi0::mainBackground =
    {"Eta3Pi0",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_3Pi0_6g)};
const std::vector<triplePi0::named_channel_t> triplePi0::otherBackgrounds =
{
    {"2Pi0Prod",    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)},
    {"EtaPi0Prod",  ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)}
};

triplePi0::triplePi0(const string& name, ant::OptionsPtr opts):
    Physics(name, opts),
    phSettings(opts->Get<bool>("AllChannels",false)),
    fitterSig("fitterSig", signal.DecayTree,         uncertModel, true ),
    fitterBkg("fitterBkg", mainBackground.DecayTree, uncertModel, true ),
    kinFitterEMB("fitterEMB", 6,                     uncertModel, true )
{
    fitterSig.SetZVertexSigma(phSettings.fitter_ZVertex);
    fitterBkg.SetZVertexSigma(phSettings.fitter_ZVertex);
    kinFitterEMB.SetZVertexSigma(phSettings.fitter_ZVertex);

    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

    promptrandom.AddPromptRange(phSettings.Range_Prompt);
    for ( const auto& range: phSettings.Ranges_Random)
        promptrandom.AddRandomRange(range);

    hist_steps      = HistFac.makeTH1D("steps","","# evts.",BinSettings(1,0,0),"steps");
    hist_channels   = HistFac.makeTH1D("channels","","# evts.",BinSettings(1,0,0),"channels");

    tree.CreateBranches(HistFac.makeTTree(phSettings.Tree_Name));
}

void triplePi0::ProcessEvent(const ant::TEvent& event, manager_t&)
{
    const auto& data   = event.Reconstructed();

    hist_steps->Fill("seen",1);

    tree.CBESum = data.Trigger.CBEnergySum;
    if (tree.CBESum < phSettings.Cut_CBESum )
        return;
    hist_steps->Fill("CB energy sum",1);


    if ( data.Candidates.size() != phSettings.Cut_NCands )
        return;
    hist_steps->Fill("N candidates",1);


    const auto& mcTrue       = event.MCTrue();
    auto& particleTree = event.MCTrue().ParticleTree;


    //===================== TreeMatching  =====================================================
    tree.MCTrue = phSettings.Index_Data;
    if (particleTree)
    {
        if (particleTree->IsEqual(signal.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            tree.MCTrue = phSettings.Index_Signal;
            hist_channels->Fill(signal.Name.c_str(),1);
        }
        else if (particleTree->IsEqual(mainBackground.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            tree.MCTrue = phSettings.Index_MainBkg;
            hist_channels->Fill(mainBackground.Name.c_str(),1);
        }
        else
        {
            auto index = phSettings.Index_Offset;
            bool found = false;
            for (const auto& otherChannel:otherBackgrounds)
            {
                if (particleTree->IsEqual(otherChannel.DecayTree,utils::ParticleTools::MatchByParticleName))
                {
                    tree.MCTrue = index;
                    hist_channels->Fill(otherChannel.Name.c_str(),1);
                    found = true;
                }
                index++;
            }
            if (!found)
            {
                tree.MCTrue = phSettings.Index_Unknown;
                hist_channels->Fill(utils::ParticleTools::GetDecayString(particleTree).c_str(),1);
            }
        }

    }




}



AUTO_REGISTER_PHYSICS(triplePi0)
