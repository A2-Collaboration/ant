#include "singlePi0.h"


#include "expconfig/ExpConfig.h"

#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/vec/LorentzVec.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const singlePi0::named_channel_t singlePi0::signal =
    {"Pi0Prod",    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)};
const singlePi0::named_channel_t singlePi0::mainBackground =
    {"EtaGG",     ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_2g)};
const std::vector<singlePi0::named_channel_t> singlePi0::otherBackgrounds =
{
    {"Pi0PiPi",      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_eeg)}
};

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
            temp+=(*(ph->Get().Leaf->Particle));
        }
        acc.push_back(temp);
    }
    return acc;
};
//auto getLorentzSumFitted = [](const vector<utils::TreeFitter::tree_t>& nodes)
//{
//    vector<TLorentzVector> acc;
//    for ( const auto& node: nodes)
//    {
//        acc.push_back(node->Get().LVSum);
//    }
//    return acc;
//};

singlePi0::singlePi0(const string& name, ant::OptionsPtr opts):
    Physics(name, opts),
    phSettings(),
    kinFitterEMB(uncertModel, true ),
    fitterSig(signal.DecayTree, uncertModel, true )
{
    fitterSig.SetZVertexSigma(phSettings.fitter_ZVertex);
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

void singlePi0::ProcessEvent(const ant::TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    const auto& data   = event.Reconstructed();

    FillStep("seen");

    tree.CBESum = triggersimu.GetCBEnergySum();

    //simulate cb-esum-trigger
    if (!triggersimu.HasTriggered())
        return;
    FillStep("Triggered");

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

    if ( data.Candidates.size() != phSettings.Cut_NCands )
        return;
    FillStep(std_ext::formatter() << "N candidates == " << phSettings.Cut_NCands);

    unique_ptr<protonSelection_t> bestSelection;

    //===================== Reconstruction ====================================================
    tree.CBAvgTime = triggersimu.GetRefTiming();

    auto bestProb_EMB  = 0.;
    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerHit));
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        tree.Tagg_Ch = static_cast<unsigned>(taggerHit.Channel);
        tree.Tagg_E  = taggerHit.PhotonEnergy;
        tree.Tagg_W  = promptrandom.FillWeight();

        auto temp_EMB_prob = -std_ext::inf;
        for ( auto i_proton: data.Candidates.get_iter())
        {
            const protonSelection_t selection(i_proton, data.Candidates,
                                              taggerHit.GetPhotonBeam(),
                                              taggerHit.PhotonEnergy);

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

            if ( EMB_result.Probability > temp_EMB_prob)
            {
                tree.SetRaw(selection);
                tree.SetEMB(kinFitterEMB,EMB_result);
                bestSelection = std_ext::make_unique<protonSelection_t>(selection);
                std_ext::copy_if_greater(bestProb_EMB,EMB_result.Probability);
            }
            FillStep(std_ext::formatter() << "EMB-Fit: prob. not in " << phSettings.Cut_EMB_Chi2);
        } // proton

        // cut on EMB
        if (!(phSettings.Cut_EMB_Chi2.Contains(tree.EMB_chi2)))
            continue;
        // fit worked at all?
        if (!bestSelection )
            continue;
        FillStep(std_ext::formatter() << "proton-selection passed");

        auto applyTreeFit = [&bestSelection](utils::TreeFitter& fitter,
                                             const std::vector<utils::TreeFitter::tree_t>& intermediates)
        {
            fitter.PrepareFits(bestSelection->Tagg_E, bestSelection->Proton, bestSelection->Photons);
            APLCON::Result_t result;
            auto best_prob = std_ext::NaN;
            fitRatings_t fr(0,0,0,{});
            while(fitter.NextFit(result))
                if (   (result.Status    == APLCON::Result_Status_t::Success)
                       && (std_ext::copy_if_greater(best_prob,result.Probability)))
                {
                    fr = fitRatings_t(best_prob,reducedChi2(result),result.NIterations,
                                      getLorentzSumUnfitted(intermediates));
                }

            return fr;
        };

        tree.SetSIG(applyTreeFit(fitterSig,pionsFitterSig));

        tree.Tree->Fill();
        hist_channels_end->Fill(trueChannel.c_str(),1);
    } // taggerHits
}

void singlePi0::ShowResult()
{
    const auto colz = drawoption("colz");

    canvas("summary") << hist_steps
                      << hist_channels
                      << hist_channels_end
                      << TTree_drawable(tree.Tree,"IM2g")
                      << TTree_drawable(tree.Tree,"EMB_chi2")
                      << samepad
                      << TTree_drawable(tree.Tree,"SIG_chi2")
                      << endc;

    auto c = canvas("SIG - cuts");
    for (auto chi2: {1.,5.,10.,15.,20.,30.,40.})
        c << TTree_drawable(tree.Tree,"MCTrue",(std_ext::formatter() << "SIG_chi2 < " << chi2));
    c << colz
      << TTree_drawable(tree.Tree,"MCTrue:SIG_chi2>>h(100,0,15,18,0,18)")
      << endc;

    canvas("pions") << colz
                    << TTree_drawable(tree.Tree,"EMB_proton.Theta() * 180/3.1415:EMB_proton.E()-938.3>>ht(100,0,1000,100,0,80)")
                    << TTree_drawable(tree.Tree,"SIG_pions.M()")
                    << endc;
}

void singlePi0::PionProdTree::SetRaw(const singlePi0::protonSelection_t& selection)
{
    proton = *selection.Proton;
    photons() = MakeTLorenz(selection.Photons);
    photonSum = selection.PhotonSum;
    IM2g = photonSum().M();
    proton_MM =   selection.PhotonBeam
                + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())
                - selection.PhotonSum;
    pg_copl   = std_ext::radian_to_degree(vec2::Phi_mpi_pi(selection.Proton->Phi()-selection.PhotonSum.Phi() - M_PI ));
    pMM_angle = std_ext::radian_to_degree(proton_MM().Angle(selection.Proton->p));
}



void singlePi0::PionProdTree::SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result)
{

    EMB_proton = *(kF.GetFittedProton());
    const auto fittedPhotons = kF.GetFittedPhotons();
    EMB_photons = MakeTLorenz(fittedPhotons);
    EMB_photonSum = accumulate(EMB_photons().begin(),EMB_photons().end(),TLorentzVector(0,0,0,0));
    EMB_IM2g = EMB_photonSum().M();
    EMB_Ebeam  = kF.GetFittedBeamE();
    EMB_iterations = result.NIterations;
    EMB_prob = result.Probability;
    EMB_chi2 = reducedChi2(result);

}

void singlePi0::PionProdTree::SetSIG(const singlePi0::fitRatings_t& fitRating)
{
    SIG_prob = fitRating.Prob;
    SIG_chi2 = fitRating.Chi2;
    SIG_iterations = fitRating.Niter;
    SIG_pions = fitRating.Intermediates;
}

AUTO_REGISTER_PHYSICS(singlePi0)
