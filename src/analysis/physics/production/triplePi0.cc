#include "triplePi0.h"

//#include "base/ParticleType.h"
//#include "base/ParticleTypeTree.h"
//#include "base/std_ext/math.h"

#include "expconfig/ExpConfig.h"

//#include "plot/root_draw.h"

#include "utils/combinatorics.h"
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
    {"2Pi04g",       ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)},
    {"EtaPi04g",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)},
    {"Eta4Pi0",      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_4Pi0_8g)},
    {"EtaPi04gPiPi", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_2gPiPi2g)},
    {"Pi0PiPi",      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0PiPi_2gPiPi)},
    {"2Pi0PiPi",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0PiPi_4gPiPi)}
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

    auto extractS = [] ( vector<utils::TreeFitter::tree_t>& nodes,
                         const utils::TreeFitter& fitter,
                         const ParticleTypeDatabase::Type& mother,
                         const ParticleTypeDatabase::Type& daughterSelection )
    {
        const auto head = fitter.GetTreeNode(mother);
        for (const auto& daughter: head->Daughters())
        {
            if (daughter->Get().TypeTree->Get() == daughterSelection)
                nodes.emplace_back(daughter);
        }
    };

    extractS(pionsFitterSig, fitterSig,
             ParticleTypeDatabase::BeamProton,
             ParticleTypeDatabase::Pi0);
    extractS(pionsFitterBkg, fitterBkg,
             ParticleTypeDatabase::Eta,
             ParticleTypeDatabase::Pi0);


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
    auto probBins = BinSettings(500,0,1);
    hist_EMB_prob       = HistFac.makeTH1D("EMB probability","Prob.","# evts.",probBins,"channels_end");
    hist_SIG_prob       = HistFac.makeTH1D("Signal-Tree probability","Prob.","# evts.",probBins,"channels_end");

    tree.CreateBranches(HistFac.makeTTree(phSettings.Tree_Name));
    tree.photons().resize(phSettings.nPhotons);
    tree.EMB_photons().resize(phSettings.nPhotons);
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


//    const auto& mcTrue       = event.MCTrue();
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

    unique_ptr<PionProdTree::particleStorage_t> bestSelection_EMB ;
    unique_ptr<PionProdTree::particleStorage_t> bestSelection_TREE;

    //===================== Reconstruction ====================================================
    tree.CBAvgTime = data.Trigger.CBTiming;

    auto bestProb_EMB = 0.;
    auto bestProb_TREE = 0.;
    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");


        promptrandom.SetTaggerHit(taggerHit.Time - tree.CBAvgTime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        tree.Tagg_Ch = static_cast<unsigned>(taggerHit.Channel);
        tree.Tagg_E  = taggerHit.PhotonEnergy;
        tree.Tagg_W  = promptrandom.FillWeight();

        tree.EMB_prob = -std_ext::inf;
        for ( auto i_proton: data.Candidates.get_iter())
        {
            auto selection = makeParticles(i_proton,data.Candidates);
            // calculate and set raw values:
            tree.SetRaw(selection,taggerHit.GetPhotonBeam());


            if (!phSettings.Cut_ProtonCopl.Contains(tree.pg_copl))
                continue;
            FillStep(std_ext::formatter() << "proton-photons coplanarity in " << phSettings.Cut_ProtonCopl);



            if ( !(phSettings.Cut_MM.Contains(tree.proton_MM().M())))
                continue;
            FillStep(std_ext::formatter() << "MM(proton) in " << phSettings.Cut_MM);

            if ( phSettings.Cut_MMAngle < (tree.pMM_angle))
                continue;
            FillStep(std_ext::formatter() << "angle(MM,proton) > " << phSettings.Cut_MMAngle);

            kinFitterEMB.SetProton(selection.Proton);
            kinFitterEMB.SetPhotons(selection.Photons);
            kinFitterEMB.SetEgammaBeam(tree.Tagg_E);

            auto EMB_result = kinFitterEMB.DoFit();
            if (!(EMB_result.Status == APLCON::Result_Status_t::Success))
                continue;

            if ( EMB_result.Probability > tree.EMB_prob)
            {
                tree.SetRaw(selection,taggerHit.GetPhotonBeam());
                tree.SetEMB(kinFitterEMB,EMB_result);
                bestSelection_EMB = std_ext::make_unique<PionProdTree::particleStorage_t>(selection);//move(selection));
                std_ext::copy_if_greater(bestProb_EMB,EMB_result.Probability);

            }
            if (!(phSettings.Cut_EMBProb.Contains(tree.EMB_prob)))
                continue;
            FillStep(std_ext::formatter() << "EMB-Fit: prob. not in " << phSettings.Cut_EMBProb);

            //==>  Signal
            fitterSig.SetProton( selection.Proton);
            fitterSig.SetPhotons(selection.Photons);
            fitterSig.SetEgammaBeam(tree.Tagg_E);

            APLCON::Result_t SIG_result;
            tree.SIG_prob = std_ext::NaN;
            while(fitterSig.NextFit(SIG_result))
            {
                if (   (SIG_result.Status    == APLCON::Result_Status_t::Success)
                       && (std_ext::copy_if_greater(tree.SIG_prob,SIG_result.Probability)))
                {
                    tree.SIG_iterations = SIG_result.NIterations;
                    tree.SetSIG(selection,taggerHit.GetPhotonBeam());
                    bestSelection_TREE = std_ext::make_unique<PionProdTree::particleStorage_t>(selection);//move(selection));
                    std_ext::copy_if_greater(bestProb_TREE,tree.SIG_prob());

                }
            }

        } // proton


        // fit worked at all?
        if (!bestSelection_EMB )
            continue;
        if (!bestSelection_TREE )
            continue;
        FillStep(std_ext::formatter() << "proton-selection passed");

        hist_EMB_prob->Fill(bestProb_EMB);
        hist_SIG_prob->Fill(bestProb_TREE);



        //==>  Background
        fitterBkg.SetProton(bestSelection_TREE->Proton);
        fitterBkg.SetPhotons(bestSelection_TREE->Photons);
        fitterBkg.SetEgammaBeam(tree.Tagg_E);

        APLCON::Result_t BKG_result;
        tree.BKG_prob = std_ext::NaN;
        while(fitterBkg.NextFit(BKG_result))
        {
            if (   (BKG_result.Status    == APLCON::Result_Status_t::Success)
                   && (std_ext::copy_if_greater(tree.BKG_prob,BKG_result.Probability)))
            {
                tree.BKG_iterations = BKG_result.NIterations;
            }
        }


        tree.Tree->Fill();
        hist_channels_end->Fill(trueChannel.c_str(),1);
    } // taggerHits

}

void triplePi0::ShowResult()
{
    canvas("summary") << hist_steps
                      << hist_channels
                      << hist_channels_end
                      << TTree_drawable(tree.Tree,"IM6g")
                      << hist_EMB_prob
                      << samepad
                      << hist_SIG_prob
                      << endc;
}



void triplePi0::PionProdTree::SetRaw(const triplePi0::PionProdTree::particleStorage_t& selection,
                                      const LorentzVec& photonBeam)
{
    proton = *selection.Proton;
    photons() = MakeTLorenz(selection.Photons);
    photonSum = selection.PhotonSum;
    IM6g = photonSum().M();
    proton_MM =   photonBeam
                + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())
                - selection.PhotonSum;
    pg_copl   = std_ext::radian_to_degree(vec2::Phi_mpi_pi(selection.Proton->Phi()-selection.PhotonSum.Phi() - M_PI ));
    pMM_angle = std_ext::radian_to_degree(proton_MM().Angle(selection.Proton->p));
}

void triplePi0::PionProdTree::SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result)
{

    EMB_proton = *(kF.GetFittedProton());
    const auto fittedPhotons = kF.GetFittedPhotons();
    EMB_photons = MakeTLorenz(fittedPhotons);
    EMB_photonSum = accumulate(EMB_photons().begin(),EMB_photons().end(),TLorentzVector(0,0,0,0));
    EMB_IM6g = EMB_photonSum().M();
    EMB_Ebeam  = kF.GetFittedBeamE();
    EMB_iterations = result.NIterations;
    EMB_prob = result.Probability;

}

void triplePi0::PionProdTree::SetSIG(const triplePi0::PionProdTree::particleStorage_t& selection, const LorentzVec& photonBeam)
{
    SIG_proton = *selection.Proton;
    SIG_photons() = MakeTLorenz(selection.Photons);
    SIG_photonSum = selection.PhotonSum;
    SIG_IM6g = SIG_photonSum().M();
    SIG_proton_MM =   photonBeam
                + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())
                - selection.PhotonSum;
    SIG_pg_copl   = std_ext::radian_to_degree(vec2::Phi_mpi_pi(selection.Proton->Phi()-selection.PhotonSum.Phi() - M_PI ));
}

AUTO_REGISTER_PHYSICS(triplePi0)
