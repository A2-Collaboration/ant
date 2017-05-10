#include "singlePi0.h"


#include "expconfig/ExpConfig.h"

#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/vec/LorentzVec.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "base/Logger.h"

#include "slowcontrol/SlowControlVariables.h"


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
    {"RhoPiPi", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Rho_PiPi)},
    {"Pi0PiPi", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_eeg)}
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
auto getLorentzSumFitted = [](const vector<utils::TreeFitter::tree_t>& nodes)
{
    vector<TLorentzVector> acc;
    for ( const auto& node: nodes)
    {
        acc.push_back(node->Get().LVSum);
    }
    return acc;
};

singlePi0::singlePi0(const string& name, ant::OptionsPtr opts):
    Physics(name, opts),
    phSettings(),
    flag_mc(opts->Get<bool>("mc", false)),
    tagger(ExpConfig::Setup::GetDetector<TaggerDetector_t>()),
    uncertModelData(// use Interpolated, based on Sergey's model
                  utils::UncertaintyModels::Interpolated::makeAndLoad(
                      utils::UncertaintyModels::Interpolated::Type_t::Data,
                      // use Sergey as starting point
                      make_shared<utils::UncertaintyModels::FitterSergey>()
                      )),
    uncertModelMC(// use Interpolated, based on Sergey's model
                utils::UncertaintyModels::Interpolated::makeAndLoad(
                    utils::UncertaintyModels::Interpolated::Type_t::MC,
                    // use Sergey as starting point
                    make_shared<utils::UncertaintyModels::FitterSergey>()
                    )),
    fitterEMB(uncertModelData, true ),
    fitterSig(signal.DecayTree, uncertModelData, true )
{
    fitterEMB.SetZVertexSigma(phSettings.fitter_ZVertex);
    fitterSig.SetZVertexSigma(phSettings.fitter_ZVertex);

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
    hist_tagger_hits    = HistFac.makeTH1D("tagger hits","# tagger hits","#",BinSettings(1,0,0),"tagger_hits");
    hist_channels       = HistFac.makeTH1D("channels","","# evts.",BinSettings(1,0,0),"channels");
    hist_channels_end   = HistFac.makeTH1D("channel-selected","","# evts.",BinSettings(1,0,0),"channels_end");

    hist_neutrals_channels
            = HistFac.makeTH2D("# neutral candidates","","# neutrals",BinSettings(1,0,0),BinSettings(5),"channels_neutrals");

    tree.CreateBranches(HistFac.makeTTree(phSettings.Tree_Name));
    tree.photons().resize(phSettings.nPhotons);
    tree.EMB_photons().resize(phSettings.nPhotons);


    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");
    const auto nchannels = Tagger->GetNChannels();


    if (!flag_mc)
    {
        slowcontrol::Variables::TaggerScalers->Request();
        slowcontrol::Variables::Trigger->Request();
    }
    fitterSig.SetUncertaintyModel(flag_mc ? uncertModelMC : uncertModelData);
    fitterEMB.SetUncertaintyModel(flag_mc ? uncertModelMC : uncertModelData);

    seenMC        = HistFac.makeTH1D("seenMC","ch","#",BinSettings(nchannels),"seenMC");
    taggerScalars = HistFac.makeTH1D("electrons","ch","# tagger scalar counts",
                                     BinSettings(nchannels),"taggerScalars",true);


    tree.TaggRates().resize(nchannels);
}

void singlePi0::ProcessEvent(const ant::TEvent& event, manager_t&)
{


    triggersimu.ProcessEvent(event);

    const auto& data   = event.Reconstructed();

    FillStep("seen");

    if ( flag_mc != data.ID.isSet(TID::Flags_t::MC))
        throw runtime_error("provided mc flag does not match input files!");

    if (flag_mc)
    {
        if (event.MCTrue().TaggerHits.size() == 1)
        {
            seenMC->Fill(event.MCTrue().TaggerHits.at(0).Channel);
        }
    }
    hist_tagger_hits->Fill(event.MCTrue().TaggerHits.size());
    if(slowcontrol::Variables::TaggerScalers->HasChanged())
    {
        const auto counts = slowcontrol::Variables::TaggerScalers->GetCounts();
        for( auto i = 0u; i < counts.size() ; ++i)
        {
            taggerScalars->Fill(i,counts.at(i));
        }
    }



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
                trueChannel = "u: " + utils::ParticleTools::GetDecayString(particleTree);
            }
        }

    }
    hist_channels->Fill(trueChannel.c_str(),1);


    tree.CBESum = triggersimu.GetCBEnergySum();

    //simulate cb-esum-trigger
    if (!triggersimu.HasTriggered())
        return;
    FillStep("Triggered");



    if (tools::cutOn("N_{cands}",phSettings.Cut_NCands,data.Candidates.size(),hist_steps))
        return;

    unique_ptr<protonSelection_t> bestSelection;

    //===================== Reconstruction ====================================================
    tree.CBAvgTime = triggersimu.GetRefTiming();

    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerHit));

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        FillStep("taggerhit inside");

        tree.Tagg_Ch = static_cast<unsigned>(taggerHit.Channel);
        tree.Tagg_E  = taggerHit.PhotonEnergy;
        tree.Tagg_W  = promptrandom.FillWeight();

        {
            const auto taggEff = tagger->GetTaggEff(taggerHit.Channel);
            tree.Tagg_Eff      = taggEff.Value;
            tree.Tagg_EffErr  = taggEff.Error;
        }



        const auto pSelections = tools::makeProtonSelections(data.Candidates,
                                                             taggerHit.GetPhotonBeam(),
                                                             taggerHit.PhotonEnergy,
                                                             phSettings.Cut_MM);

        auto temp_prob = 0.;
        auto bestFound = false;
        for ( const auto& selection: pSelections)
        {

            const auto EMB_result = fitterEMB.DoFit(selection.Tagg_E, selection.Proton, selection.Photons);
            if (EMB_result.Status != APLCON::Result_Status_t::Success)
                continue;
            if (tools::cutOn("EMB-prob && success",phSettings.Cut_EMB_prob,EMB_result.Probability,hist_steps))
                continue;
            const auto sigFitRatings
                    = [&selection](utils::TreeFitter& fitter,
                                   const std::vector<utils::TreeFitter::tree_t>& intermediates)
                      {
                          fitter.PrepareFits(selection.Tagg_E, selection.Proton, selection.Photons);
                          APLCON::Result_t result;
                          auto best_prob = std_ext::NaN;
                          fitRatings_t fr(0,0,0,false,{});
                          while(fitter.NextFit(result))
                              if (   (result.Status    == APLCON::Result_Status_t::Success)
                                     && (std_ext::copy_if_greater(best_prob,result.Probability)))
                              {
                                  fr = fitRatings_t(best_prob,reducedChi2(result),result.NIterations,
                                                    result.Status == APLCON::Result_Status_t::Success,
                                                    getLorentzSumFitted(intermediates));
                              }

                          return fr;
                      }(fitterSig,pionsFitterSig);

            if (!(sigFitRatings.FitOk)) continue;
            FillStep(std_ext::formatter() << "Fit succesful");



            if ( sigFitRatings.Prob > temp_prob)
            {
                bestFound = true;
                temp_prob = sigFitRatings.Prob;
                tree.SetRaw(selection);
                tree.SetEMB(fitterEMB,EMB_result);
                tree.SetSIG(sigFitRatings);
                tree.ProtonVetoE() = selection.ProtonVetoE;
                tree.PionVetoE()   = selection.PhotonVetoE;
            }

        } // proton

        if (bestFound)
        {
            const auto eThresh = 0.0;
            const auto neutral_cands = tools::getNeutral(data,eThresh);

            tree.Neutrals() = neutral_cands.size();

            tree.Tree->Fill();
            hist_channels_end->Fill(trueChannel.c_str(),1);
            hist_neutrals_channels->Fill(trueChannel.c_str(),neutral_cands.size(),1);

            if (flag_mc)
            {
                for (auto& tr: tree.TaggRates())
                    tr = 1;
                tree.ExpLivetime() = 1;
            }
            else if(slowcontrol::Variables::TaggerScalers->HasChanged())
            {
                tree.TaggRates()  = slowcontrol::Variables::TaggerScalers->Get();
                tree.ExpLivetime()  = slowcontrol::Variables::Trigger->GetExpLivetime();
            }
        }
    } // taggerHits
}

void singlePi0::ShowResult()
{
    const auto colz = drawoption("colz");

    canvas("summary") << hist_steps
                      << hist_channels
                      << hist_channels_end
                      << endc;

    canvas("channels") << colz
            << hist_neutrals_channels
            << endc;
    canvas("norm") << seenMC << taggerScalars << endc;
}

void singlePi0::PionProdTree::SetRaw(const tools::protonSelection_t& selection)
{
    proton = *selection.Proton;
    protonTime = selection.Proton->Candidate->Time;

    photons() = tools::MakeTLorenz(selection.Photons);
    photonTimes().resize(selection.Photons.size());
    transform(selection.Photons.begin(),selection.Photons.end(),
              photonTimes().begin(),
              [](const TParticlePtr& photon) {return photon->Candidate->Time;});


    photonSum = selection.PhotonSum;
    IM2g      = photonSum().M();
    proton_MM = selection.Proton_MM;
    pg_copl    = selection.Copl_pg;
    pMM_angle  = selection.Angle_pMM;
}



void singlePi0::PionProdTree::SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result)
{
    const auto fittedPhotons = kF.GetFittedPhotons();

    EMB_proton     = *(kF.GetFittedProton());
    EMB_photons    = tools::MakeTLorenz(fittedPhotons);
    EMB_photonSum  = accumulate(EMB_photons().begin(),EMB_photons().end(),TLorentzVector(0,0,0,0));
    EMB_IM2g       = EMB_photonSum().M();
    EMB_Ebeam      = kF.GetFittedBeamE();
    EMB_iterations = result.NIterations;
    EMB_prob       = result.Probability;
    EMB_chi2       = reducedChi2(result);

}

void singlePi0::PionProdTree::SetSIG(const singlePi0::fitRatings_t& fitRating)
{
    SIG_prob = fitRating.Prob;
    SIG_chi2 = fitRating.Chi2;
    SIG_iterations = fitRating.Niter;
    SIG_pions = fitRating.Intermediates;
}


AUTO_REGISTER_PHYSICS(singlePi0)
