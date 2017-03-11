#include "triplePi0.h"

//#include "base/ParticleType.h"
//#include "base/ParticleTypeTree.h"
//#include "base/std_ext/math.h"

#include "expconfig/ExpConfig.h"

//#include "plot/root_draw.h"

#include "utils/combinatorics.h"
#include "utils/particle_tools.h"
#include "base/vec/LorentzVec.h"
#include "base/Logger.h"
//#include "utils/ParticleID.h"

//#include <algorithm>
//#include <cassert>
//#include <chrono>


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

#include "utils/uncertainties/Interpolated.h"


const triplePi0::named_channel_t triplePi0::signal =
    {"3Pi0Prod",    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)};
const triplePi0::named_channel_t triplePi0::mainBackground =
    {"Eta3Pi0",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_3Pi0_6g)};
const triplePi0::named_channel_t triplePi0::sigmaBackground =
    {"SigmaK0S",   ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::SigmaPlusK0s_6g)};
const std::vector<triplePi0::named_channel_t> triplePi0::otherBackgrounds =
{
    {"2Pi04g",       ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)},
    {"EtaPi04g",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)},
    {"Eta4Pi0",      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_4Pi0_8g)},
    {"EtaPi04gPiPi", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_2gPiPi2g)},
    {"Pi0PiPi",      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0PiPi_2gPiPi)},
    {"2Pi0PiPi",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0PiPi_4gPiPi)},
    {"Etap3Pi0",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_3Pi0_6g)}
};

string triplePi0::getOtherChannelNames(const unsigned i)
{
    if (i == 9)
        return "unkown";
    if (i>=10 && i<10+otherBackgrounds.size())
            return otherBackgrounds.at(i-10).Name;
    return "error";
}

auto reducedChi2 = [](const APLCON::Result_t& ares)
{
    return 1. * ares.ChiSquare / ares.NDoF;
};

auto getLorentzSumUnfitted = [](const vector<utils::TreeFitter::tree_t>& nodes)
{
    vector<TLorentzVector> acc;
    for ( const auto& node: nodes)
    {
        LorentzVec temp({0,0,0},0);
        for ( const auto& ph: node->Daughters())
        {
            temp+=(*(ph->Get().Leave->Particle));
        }
        acc.push_back(temp);
    }
    return acc;
};

auto getTLorentz = [] (const utils::TreeFitter::tree_t& node)
{
    return node->Get().LVSum;
};

auto getTreeFitPhotonIndices = [] (const TParticleList& orig_Photons,
                                   const utils::TreeFitter& treeFitter)
{
    const auto allP = treeFitter.GetFitParticles();
    const vector<utils::Fitter::FitParticle> fitPhotons(allP.begin()+1,
                                                        allP.end());

    vector<unsigned> combination;
    for (unsigned iFit = 0; iFit < fitPhotons.size(); ++iFit)
    {
        for (unsigned iOrig = 0; iOrig < orig_Photons.size(); ++iOrig)
        {
            if ( fitPhotons.at(iFit).Particle == orig_Photons.at(iOrig))
            {
                combination.push_back(iOrig);
                continue;
            }
        }
    }
    return combination;
};

auto getLorentzSumFitted = [](const vector<utils::TreeFitter::tree_t>& nodes)
{
    vector<TLorentzVector> acc;
    for ( const auto& node: nodes)
    {
        acc.push_back(node->Get().LVSum);
    }
    return acc;
};



triplePi0::triplePi0(const string& name, ant::OptionsPtr opts):
    Physics(name, opts),
    phSettings(),
    tagger(ExpConfig::Setup::GetDetector<TaggerDetector_t>()),
    uncertModel(utils::UncertaintyModels::Interpolated::makeAndLoad()),
    kinFitterEMB("fitterEMB", 6,                                  uncertModel, true ),
    fitterSig("fitterSig", signal.DecayTree,                      uncertModel, true ),
    fitterBkg("fitterBkg", mainBackground.DecayTree,              uncertModel, true ),
    fitterSigmaPlus("fittedSigmaPlus", sigmaBackground.DecayTree, uncertModel, true )
{
    fitterSig.SetZVertexSigma(phSettings.fitter_ZVertex);
    fitterBkg.SetZVertexSigma(phSettings.fitter_ZVertex);
    fitterSigmaPlus.SetZVertexSigma(phSettings.fitter_ZVertex);
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
    pionsFitterSigmaPlus = fitterSigmaPlus.GetTreeNodes(ParticleTypeDatabase::Pi0);

    kaonFitterSigmaPlus  = fitterSigmaPlus.GetTreeNode(ParticleTypeDatabase::K0s);
    sigmaFitterSigmaPlus = fitterSigmaPlus.GetTreeNode(ParticleTypeDatabase::SigmaPlus);

    // be lazy and catch complete class...
    fitterSigmaPlus.SetIterationFilter([this] () {
        const auto sigmaPlus_cut = ParticleTypeDatabase::SigmaPlus.GetWindow(200);
        const auto K0s_cut = ParticleTypeDatabase::K0s.GetWindow(100);
        auto ok = sigmaPlus_cut.Contains(sigmaFitterSigmaPlus->Get().LVSum.M()) &&
                  K0s_cut.Contains(kaonFitterSigmaPlus->Get().LVSum.M());
        return ok;
    });


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
    tree.photons().resize(phSettings.nPhotons);
    tree.EMB_photons().resize(phSettings.nPhotons);
}

const triplePi0::fitRatings_t applyTreeFit(utils::TreeFitter& fitter,
                                           const std::vector<utils::TreeFitter::tree_t>& intermediates,
                                           const triplePi0::protonSelection_t& protonSelection)
{
    fitter.PrepareFits(protonSelection.Tagg_E,
                       protonSelection.Proton,
                       protonSelection.Photons);
    APLCON::Result_t result;
    auto best_prob = std_ext::NaN;
    triplePi0::fitRatings_t fr(0,0,0,
                               {0,0,0,0},
                               {},{});
    while(fitter.NextFit(result))
        if (   (result.Status    == APLCON::Result_Status_t::Success)
               && (std_ext::copy_if_greater(best_prob,result.Probability)))
        {
            fr = triplePi0::fitRatings_t(best_prob,reducedChi2(result),result.NIterations,
                                         *fitter.GetFittedProton(),
                                         getLorentzSumFitted(intermediates),
                                         getTreeFitPhotonIndices(protonSelection.Photons,fitter));
        }

    return fr;
}

void triplePi0::ProcessEvent(const ant::TEvent& event, manager_t&)
{
    const auto& data   = event.Reconstructed();

    FillStep("seen");

    tree.CBESum = data.Trigger.CBEnergySum;

    //simulate cb-esum-trigger
    if (tree.CBESum < phSettings.Cut_CBESum )
        return;
    FillStep(std_ext::formatter() << "CB energy sum < " << phSettings.Cut_CBESum);

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
        else if (particleTree->IsEqual(sigmaBackground.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            tree.MCTrue = phSettings.Index_SigmaBkg;
            trueChannel = sigmaBackground.Name;
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
                trueChannel = utils::ParticleTools::GetDecayString(particleTree) + ": unknown";
            }
        }

    }
    hist_channels->Fill(trueChannel.c_str(),1);

    if ( data.Candidates.size() != phSettings.Cut_NCands )
        return;
    FillStep(std_ext::formatter() << "N candidates == " << phSettings.Cut_NCands);


    //===================== Reconstruction ====================================================
    tree.CBAvgTime = data.Trigger.CBTiming;

    auto bestProb_SIG  = 0.;
    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");

        promptrandom.SetTaggerHit(taggerHit.Time - tree.CBAvgTime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        tree.Tagg_Ch  = static_cast<unsigned>(taggerHit.Channel);
        tree.Tagg_E   = taggerHit.PhotonEnergy;
        tree.Tagg_W   = promptrandom.FillWeight();

        {
            const auto taggEff = tagger->GetTaggEff(taggerHit.Channel);
            tree.Tagg_Eff      = taggEff.Value;
            tree.Tagg_EffErr  = taggEff.Error;
        }

        for ( auto i_proton: data.Candidates.get_iter())
        {
            const protonSelection_t selection(i_proton, data.Candidates,
                                              taggerHit.GetPhotonBeam(),
                                              taggerHit.PhotonEnergy);




            // cuts "to save CPU time"
            if (!phSettings.Cut_ProtonCopl.Contains(selection.Copl_pg))
                continue;
            FillStep(std_ext::formatter() << "proton-photons coplanarity in " << phSettings.Cut_ProtonCopl);

            if ( !(phSettings.Cut_MM.Contains(selection.Proton_MM.M())))
                continue;
            FillStep(std_ext::formatter() << "MM(proton) in " << phSettings.Cut_MM);

            if ( phSettings.Cut_MMAngle < (selection.Angle_pMM))
                continue;
            FillStep(std_ext::formatter() << "angle(MM,proton) > " << phSettings.Cut_MMAngle);


            auto EMB_result = kinFitterEMB.DoFit(selection.Tagg_E, selection.Proton, selection.Photons);
            if (!(EMB_result.Status == APLCON::Result_Status_t::Success))
                continue;
            FillStep(std_ext::formatter() << "EMB-prefit succesful");

            if ( EMB_result.Probability < phSettings.Cut_EMB_prob)
                continue;
            FillStep(std_ext::formatter() << "EMB-prob > " << phSettings.Cut_EMB_prob);


            // let signal-tree-fitter decide about the right comination
            auto sigFitRatings = applyTreeFit(fitterSig,pionsFitterSig,selection);
            if ( sigFitRatings.Prob > bestProb_SIG )
            {
                bestProb_SIG = sigFitRatings.Prob;

                tree.SetRaw(selection);
                tree.SetEMB(kinFitterEMB,EMB_result);
                tree.SetSIG(sigFitRatings);
                tree.SetBKG(applyTreeFit(fitterBkg,pionsFitterBkg,selection));
                // do sigmaplus-k0 treefit
                /*
                {
                    APLCON::Result_t result;
                    auto best_prob = std_ext::NaN;


                    fitterSigmaPlus.PrepareFits(selection.Tagg_E,
                                                selection.Proton,
                                                selection.Photons);

                    while(fitterSigmaPlus.NextFit(result))
                    {
                        if ( (result.Status == APLCON::Result_Status_t::Success)
                             && (std_ext::copy_if_greater(best_prob,result.Probability)))
                        {
                            tree.SIGMA_prob = best_prob;
                            tree.SIGMA_chi2 = reducedChi2(result);
                            tree.SIGMA_combination() = (getTreeFitPhotonIndices(selection.Photons,
                                                                                fitterSigmaPlus));
                            tree.SIGMA_pions()   = getLorentzSumFitted(pionsFitterSigmaPlus);
                            tree.SIGMA_k0s       = getTLorentz(kaonFitterSigmaPlus);
                            tree.SIGMA_SigmaPlus = getTLorentz(sigmaFitterSigmaPlus);
                        }
                    }
                }  // scope for SIGMA - treefit
                */
            } // endif best SIG - treefit

        } // proton - candidate - loop

        tree.Tree->Fill();
        hist_channels_end->Fill(trueChannel.c_str(),1);

    } // taggerHits - loop
}

void triplePi0::ShowResult()
{

    canvas("summary") << hist_steps
                      << hist_channels
                      << hist_channels_end
                      << TTree_drawable(tree.Tree,"IM6g")
                      << endc;
}

void triplePi0::PionProdTree::SetRaw(const triplePi0::protonSelection_t& selection)
{
    proton     = *selection.Proton;
    protonTime = selection.Proton->Candidate->Time;

    photons()  = MakeTLorenz(selection.Photons);
    photonTimes().resize(selection.Photons.size());
    transform(selection.Photons.begin(),selection.Photons.end(),
              photonTimes().begin(),
              [](const TParticlePtr& photon) {return photon->Candidate->Time;});

    photonSum  = selection.PhotonSum;
    IM6g       = photonSum().M();
    proton_MM  = selection.Proton_MM;
    pg_copl    = selection.Copl_pg;
    pMM_angle  = selection.Angle_pMM;
}



void triplePi0::PionProdTree::SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result)
{

    const auto fittedPhotons = kF.GetFittedPhotons();
    const auto phE           = kF.GetFittedBeamE();

    EMB_proton     = *(kF.GetFittedProton());
    EMB_photons    = MakeTLorenz(fittedPhotons);
    EMB_photonSum  = accumulate(EMB_photons().begin(),EMB_photons().end(),LorentzVec({0,0,0},0));
    EMB_IM6g       = EMB_photonSum().M();
    EMB_Ebeam      = phE;

    EMB_proton_MM  =   LorentzVec({0,0,phE},phE) + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())
                     - EMB_photonSum();

    EMB_iterations = result.NIterations;
    EMB_prob       = result.Probability;
    EMB_chi2       = reducedChi2(result);

}

void triplePi0::PionProdTree::SetSIG(const triplePi0::fitRatings_t& fitRating)
{
    SIG_prob        = fitRating.Prob;
    SIG_chi2        = fitRating.Chi2;
    SIG_iterations  = fitRating.Niter;
    SIG_pions       = fitRating.Intermediates;
    SIG_IM3Pi0      = accumulate(fitRating.Intermediates.begin(),
                                 fitRating.Intermediates.end(),
                                 LorentzVec({0,0,0},0)).M();
    SIG_proton      = fitRating.Proton;
    SIG_combination = fitRating.PhotonCombination;
}

void triplePi0::PionProdTree::SetBKG(const triplePi0::fitRatings_t& fitRating)
{
    BKG_prob        = fitRating.Prob;
    BKG_chi2        = fitRating.Chi2;
    BKG_iterations  = fitRating.Niter;
    BKG_pions       = fitRating.Intermediates;
    BKG_combination = fitRating.PhotonCombination;
}




AUTO_REGISTER_PHYSICS(triplePi0)
