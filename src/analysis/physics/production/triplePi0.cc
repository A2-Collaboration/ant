#include "triplePi0.h"

//#include "base/ParticleType.h"
//#include "base/ParticleTypeTree.h"
//#include "base/std_ext/math.h"

#include "expconfig/ExpConfig.h"

//#include "plot/root_draw.h"

//#include "utils/combinatorics.h"
#include "utils/particle_tools.h"
#include "base/vec/LorentzVec.h"
//#include "utils/ParticleID.h"

//#include <algorithm>
//#include <cassert>
//#include <chrono>


using namespace std;
using namespace ant;
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
    kinFitterEMB("fitterEMB", 6,                     uncertModel, true ),
    fitterSig("fitterSig", signal.DecayTree,         uncertModel, true ),
    fitterBkg("fitterBkg", mainBackground.DecayTree, uncertModel, true )
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

    hist_steps          = HistFac.makeTH1D("steps","","# evts.",BinSettings(1,0,0),"steps");
    hist_channels       = HistFac.makeTH1D("channels","","# evts.",BinSettings(1,0,0),"channels");
    hist_channels_end   = HistFac.makeTH1D("channel-selected","","# evts.",BinSettings(1,0,0),"channels_end");

    tree.CreateBranches(HistFac.makeTTree(phSettings.Tree_Name));
}

void triplePi0::ProcessEvent(const ant::TEvent& event, manager_t&)
{
    const auto& data   = event.Reconstructed();

    FillStep("seen");

    tree.CBESum = data.Trigger.CBEnergySum;
    if (tree.CBESum < phSettings.Cut_CBESum )
        return;
    FillStep(std_ext::formatter() << "CB energy sum < " << phSettings.Cut_CBESum);


    if ( data.Candidates.size() != phSettings.Cut_NCands )
        return;
    FillStep(std_ext::formatter() << "N candidates == " << phSettings.Cut_NCands);


    const auto& mcTrue       = event.MCTrue();
    auto& particleTree = event.MCTrue().ParticleTree;


    //===================== TreeMatching   ====================================================
    tree.MCTrue = phSettings.Index_Data;
    string trueChannel = "Unknown/Data";
    if (particleTree)
    {
        if (particleTree->IsEqual(signal.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            tree.MCTrue = phSettings.Index_Signal;
            trueChannel = signal.Name;
        }
        else if (particleTree->IsEqual(mainBackground.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            tree.MCTrue = phSettings.Index_MainBkg;
            trueChannel = mainBackground.Name;
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
                    trueChannel = otherChannel.Name;
                    found = true;
                }
                index++;
            }
            if (!found)
            {
                tree.MCTrue = phSettings.Index_Unknown;
                trueChannel = utils::ParticleTools::GetDecayString(particleTree);
            }
        }

    }
    hist_channels->Fill(trueChannel.c_str(),1);

    //===================== Reconstruction ====================================================
    tree.CBAvgTime = data.Trigger.CBTiming;
    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");

        promptrandom.SetTaggerHit(taggerHit.Time - tree.CBAvgTime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        tree.Tagg_Ch = static_cast<unsigned>(taggerHit.Channel);
        tree.Tagg_E  = taggerHit.PhotonEnergy;
        tree.Tagg_W  = promptrandom.FillWeight();

        for ( auto i_proton: data.Candidates.get_iter())
        {
            const auto proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, i_proton);
            vector<TParticlePtr> photons;
            for ( auto i_photon : data.Candidates.get_iter())
                if (!(i_photon == i_proton))
                    photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));


            LorentzVec photonSum({0,0,0},0);
            for(const auto& p : photons) {
                photonSum += *p;
            }

            const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photonSum.Phi() - M_PI ));
            if (!phSettings.Cut_ProtonCopl.Contains(d_phi))
                continue;
            FillStep(std_ext::formatter() << "proton-photons coplanarity in " << phSettings.Cut_ProtonCopl);


            const LorentzVec missing = taggerHit.GetPhotonBeam()
                                       + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())
                                       - photonSum;
            if ( !(phSettings.Cut_MM.Contains(missing.M())))
                continue;
            FillStep(std_ext::formatter() << "MM(proton) in " << phSettings.Cut_MM);

            if ( phSettings.Cut_MMAngle < (std_ext::radian_to_degree(missing.Angle(proton->p)) ))
                continue;
            FillStep(std_ext::formatter() << "angle(MM,proton) > " << phSettings.Cut_MMAngle);


    //===================== EMB-Fitting    ====================================================
            kinFitterEMB.SetProton(proton);
            kinFitterEMB.SetPhotons(photons);
            kinFitterEMB.SetEgammaBeam(tree.Tagg_E);

            auto EMB_result = kinFitterEMB.DoFit();
            if (!(EMB_result.Status == APLCON::Result_Status_t::Success))
                continue;
            FillStep("EMB-Fit: success");
            if (!(phSettings.Cut_EMBProb.Contains(EMB_result.Probability)))
                continue;
            FillStep(std_ext::formatter() << "EMB-Fit: prob. not in " << phSettings.Cut_EMBProb);

  //===================== tree-fitting   ====================================================
            //==>  Signal

            fitterSig.SetProton(proton);
            fitterSig.SetPhotons(photons);
            fitterSig.SetEgammaBeam(tree.Tagg_E);

            APLCON::Result_t SIG_result;
            tree.SIG_prob = std_ext::NaN;
            while(fitterSig.NextFit(SIG_result))
            {
                if (   (SIG_result.Status    == APLCON::Result_Status_t::Success)
                    && (std_ext::copy_if_greater(tree.SIG_prob,SIG_result.Probability)))
                {
                    //do stuff
                }
            }

            //==>  Background

            fitterBkg.SetProton(proton);
            fitterBkg.SetPhotons(photons);
            fitterBkg.SetEgammaBeam(tree.Tagg_E);

            APLCON::Result_t BKG_result;
            tree.BKG_prob = std_ext::NaN;
            while(fitterSig.NextFit(BKG_result))
            {
                if (   (BKG_result.Status    == APLCON::Result_Status_t::Success)
                    && (std_ext::copy_if_greater(tree.BKG_prob,BKG_result.Probability)))
                {
                    //do stuff
                }
            }
        }// proton-candidates

        tree.Tree->Fill();
        hist_channels_end->Fill(trueChannel.c_str(),1);
    } // taggerHits

}


void triplePi0::ShowResult()
{
    canvas("summary") << hist_steps
                      << hist_channels
                      << hist_channels_end
                      << endc;
}

AUTO_REGISTER_PHYSICS(triplePi0)
