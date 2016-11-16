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

    tree.CreateBranches(HistFac.makeTTree(phSettings.Tree_Name));
    tree.photons().resize(phSettings.nPhotons);
    tree.EMB_photons().resize(phSettings.nPhotons);
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
    tree.CBAvgTime = data.Trigger.CBTiming;

    auto bestProb_EMB  = 0.;
    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");

        promptrandom.SetTaggerHit(taggerHit.Time - tree.CBAvgTime);
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

            kinFitterEMB.SetProton(selection.Proton);
            kinFitterEMB.SetPhotons(selection.Photons);
            kinFitterEMB.SetEgammaBeam(selection.Tagg_E);

            auto EMB_result = kinFitterEMB.DoFit();
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
            fitter.SetProton( bestSelection->Proton);
            fitter.SetPhotons(bestSelection->Photons);
            fitter.SetEgammaBeam(bestSelection->Tagg_E);
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
        tree.SetBKG(applyTreeFit(fitterBkg,pionsFitterBkg));

        tree.Tree->Fill();
        hist_channels_end->Fill(trueChannel.c_str(),1);
    } // taggerHits
}

void triplePi0::ShowResult()
{
    const auto colz = drawoption("colz");

    canvas("summary") << hist_steps
                      << hist_channels
                      << hist_channels_end
                      << TTree_drawable(tree.Tree,"IM6g")
                      << TTree_drawable(tree.Tree,"EMB_chi2")
                      << samepad
                      << TTree_drawable(tree.Tree,"SIG_chi2")
                      << samepad
                      << TTree_drawable(tree.Tree,"BKG_chi2")
                      << endc;

    auto c = canvas("SIG - cuts");
    for (auto chi2: {1.,5.,10.,15.,20.,30.,40.})
        c << TTree_drawable(tree.Tree,"MCTrue",(std_ext::formatter() << "SIG_chi2 < " << chi2));
    c << colz
      << TTree_drawable(tree.Tree,"MCTrue:SIG_chi2>>h(100,0,15,18,0,18)")
      << endc;

    canvas("SIG-BKG") << colz
                      << TTree_drawable(tree.Tree,"BKG_chi2:SIG_chi2","SIG_chi2 < 40")
                      << endc;
    canvas("pions") << colz
                    << TTree_drawable(tree.Tree,"SIG_pions.M()")
                    << samepad
                    << TTree_drawable(tree.Tree,"BKG_pions.M()")
                    << endc;
}

void triplePi0::PionProdTree::SetRaw(const triplePi0::protonSelection_t& selection)
{
    proton = *selection.Proton;
    photons() = MakeTLorenz(selection.Photons);
    photonSum = selection.PhotonSum;
    IM6g = photonSum().M();
    proton_MM =   selection.PhotonBeam
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
    EMB_chi2 = reducedChi2(result);

}

void triplePi0::PionProdTree::SetSIG(const triplePi0::fitRatings_t& fitRating)
{
    SIG_prob = fitRating.Prob;
    SIG_chi2 = fitRating.Chi2;
    SIG_iterations = fitRating.Niter;
    SIG_pions = fitRating.Intermediates;
}

void triplePi0::PionProdTree::SetBKG(const triplePi0::fitRatings_t& fitRating)
{
    BKG_prob = fitRating.Prob;
    BKG_chi2 = fitRating.Chi2;
    BKG_iterations = fitRating.Niter;
    BKG_pions = fitRating.Intermediates;
}


AUTO_REGISTER_PHYSICS(triplePi0)
