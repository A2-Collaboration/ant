#include "singlePi0.h"


#include "expconfig/ExpConfig.h"

#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/vec/LorentzVec.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "utils/ParticleTools.h"

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


auto getPi0COMS = [] (const double taggE, const LorentzVec& pi0)
{
    TLorentzVector bt(0,0, taggE, taggE + ParticleTypeDatabase::Proton.Mass());
    auto comsPi0 = pi0;
    comsPi0.Boost(-bt.BoostVector());
    return comsPi0;
};

auto getTruePi0 = [] (const TParticleTree_t& tree)
{
    return LorentzVec(*utils::ParticleTools::FindParticle(ParticleTypeDatabase::Pi0,tree));
};

auto getEgamma = [] (const TParticleTree_t& tree)
{
    return tree->Get()->Ek();
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
    fitterEMB(uncertModelData, true)
{
    fitterEMB.SetZVertexSigma(phSettings.fitter_ZVertex);


    promptrandom.AddPromptRange(phSettings.Range_Prompt);
    for ( const auto& range: phSettings.Ranges_Random)
        promptrandom.AddRandomRange(range);

    hist_steps          = HistFac.makeTH1D("steps","","# evts.",BinSettings(1,0,0),"steps");
    hist_tagger_hits    = HistFac.makeTH1D("tagger hits","# tagger hits","#",BinSettings(1,0,0),"tagger_hits");
    hist_channels       = HistFac.makeTH1D("channels","","# evts.",BinSettings(1,0,0),"channels");
    hist_channels_end   = HistFac.makeTH1D("channel-selected","","# evts.",BinSettings(1,0,0),"channels_end");
    hist_ncands         = HistFac.makeTH1D("","# candidates","# evts.",BinSettings(7),"ncands");

    hist_neutrals_channels
            = HistFac.makeTH2D("# neutral candidates","","# neutrals",BinSettings(1,0,0),BinSettings(5),"channels_neutrals");

    tree.CreateBranches(HistFac.makeTTree(phSettings.Tree_Name));
    tree.photons().resize(phSettings.nPhotons);
    tree.EMB_photons().resize(phSettings.nPhotons);

    seenSignal.CreateBranches(HistFac.makeTTree(seenSignal.treeName()));
    recSignal.CreateBranches(HistFac.makeTTree(recSignal.treeName()));

    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");
    const auto nchannels = Tagger->GetNChannels();


    if (!flag_mc)
    {
        slowcontrol::Variables::TaggerScalers->Request();
        slowcontrol::Variables::Trigger->Request();
    }
    fitterEMB.SetUncertaintyModel(flag_mc ? uncertModelMC : uncertModelData);

    seenMC        = HistFac.makeTH1D("seenMC","ch","#",BinSettings(nchannels),"seenMC");
    taggerScalars = HistFac.makeTH1D("electrons","ch","# tagger scalar counts",
                                     BinSettings(nchannels),"taggerScalars",true);


    const BinSettings taggerbins(nchannels);
    const BinSettings cosTheta(32,-1,1);

    hist_seen          = HistFac.makeTH2D("seen #pi^{0}", "taggerCh","cos(#theta_{#pi^{0}})",
                                          taggerbins, cosTheta,
                                          "seen2d",true);
    hist_rec           = HistFac.makeTH2D("reconstructed #pi^{0}", "taggerCh","cos(#theta_{#pi^{0}})",
                                          taggerbins, cosTheta,
                                          "rec2d",true);
    hist_efficiency    = HistFac.makeTH2D("efficiency", "taggerCh","cos(#theta_{#pi^{0}})",
                                          taggerbins, cosTheta,
                                          "eff2d",true);
}

void singlePi0::ProcessEvent(const ant::TEvent& event, manager_t&)
{
    const auto& data   = event.Reconstructed();

    // simulate trigger:
    triggersimu.ProcessEvent(event);
    // check if mc-flag ist set properly:
    if ( flag_mc != data.ID.isSet(TID::Flags_t::MC))
        throw runtime_error("provided mc flag does not match input files!");

    FillStep("seen");
    hist_ncands->Fill(data.Candidates.size());
    hist_tagger_hits->Fill(event.MCTrue().TaggerHits.size());


    //===================== TreeMatching   ==================================================================

    string trueChannel = "data";
    tree.MCTrue = phSettings.Index_Data;
    if (flag_mc)
    {
        tree.MCTrue() = phSettings.Index_brokenTree;
        trueChannel = "no Pluto tree";
    }
    auto& particleTree = event.MCTrue().ParticleTree;
    if (particleTree)
    {
        if (particleTree->IsEqual(signal.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            const auto& taggerhits = event.MCTrue().TaggerHits;
            const auto  truePi0    = getTruePi0(particleTree);

            if ( taggerhits.size() > 1)
            {
                LOG(INFO) << event ;
                throw runtime_error("mc should always have no more than one Taggerbin, check mctrue!");
            }

            // pluto reader generates only taggerhits, if beamtarget->Ek() is within
            // [Etagg_min,Etagg_max]!
            // see PlutoReader::CopyPluto(TEventData& mctrue)!!!
            if ( taggerhits.size() == 1)
            {
                seenSignal.Egamma()      = getEgamma(particleTree);
                seenSignal.TaggerBin()   = taggerhits[0].Channel;
                seenSignal.CosThetaPi0() = cos(getPi0COMS(seenSignal.Egamma(), truePi0).Theta());
                seenSignal.Phi()         = truePi0.Phi();
                seenSignal.Theta()       = truePi0.Theta();
                seenSignal.Tree->Fill();
            }
            tree.MCTrue() = phSettings.Index_Signal;
            trueChannel = signal.Name;
        }
        else if (particleTree->IsEqual(mainBackground.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            tree.MCTrue() = phSettings.Index_MainBkg;
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
                    tree.MCTrue() = index;
                    trueChannel = otherChannel.Name;
                    found = true;
                }
                index++;
            }
            if (!found)
            {
                tree.MCTrue() = phSettings.Index_UnTagged;
                trueChannel = utils::ParticleTools::GetDecayString(particleTree);
            }
        }

    }
    hist_channels->Fill(trueChannel.c_str(),1);



    tree.CBESum = triggersimu.GetCBEnergySum();

    //simulate cb-esum-trigger
    if (!triggersimu.HasTriggered())
        return;
    FillStep("Triggered");


    // cut on n candidates
    if (tools::cutOn("N_{cands}",phSettings.Cut_NCands,data.Candidates.size(),hist_steps))
        return;

    //===================== Reconstruction ====================================================
    tree.CBAvgTime = triggersimu.GetRefTiming();
    utils::ProtonPhotonCombs proton_photons(data.Candidates);

    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerHit));

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        FillStep("taggerhit inside");

        tree.Tagg_Ch() = static_cast<unsigned>(taggerHit.Channel);
        tree.Tagg_E()  = taggerHit.PhotonEnergy;
        tree.Tagg_W()  = promptrandom.FillWeight();

        {
            const auto taggEff = tagger->GetTaggEff(taggerHit.Channel);
            tree.Tagg_Eff()      = taggEff.Value;
            tree.Tagg_EffErr()   = taggEff.Error;
        }

//        const auto pSelections = tools::makeProtonSelections(data.Candidates,
//                                                             taggerHit.GetPhotonBeam(),
//                                                             taggerHit.PhotonEnergy,
//                                                             phSettings.Cut_MM);
        auto selections =  proton_photons()
                           .FilterMult(phSettings.nPhotons,100)
                           .FilterIM(phSettings.Cut_IM)
                           .FilterMM(taggerHit, phSettings.Cut_MM);

        if (selections.empty())
        {
            FillStep("No combs left");
            continue;
        }


        auto temp_prob = 0.;
        auto bestFound = false;
        for ( const auto& selection: selections)
        {

            const auto EMB_result = fitterEMB.DoFit(taggerHit.PhotonEnergy, selection.Proton, selection.Photons);
            if (EMB_result.Status != APLCON::Result_Status_t::Success)
                continue;

            FillStep("kin fit ok");


            if ( EMB_result.Probability > temp_prob)
            {
                bestFound = true;
                temp_prob = EMB_result.Probability;
                tree.SetRaw(selection);
                tree.SetEMB(fitterEMB,EMB_result);
                tree.PionVetoE()   = selection.Proton->Candidate->VetoEnergy;
            }

        } // proton

        if (bestFound)
        {
            FillStep("p identified");
            const auto eThresh = 0.0;
            const auto neutral_cands = tools::getNeutral(data,eThresh);

            tree.Neutrals() = neutral_cands.size();

            hist_channels_end->Fill(trueChannel.c_str(),1);
            hist_neutrals_channels->Fill(trueChannel.c_str(),neutral_cands.size(),1);

            if (flag_mc)
            {
                tree.ExpLivetime() = 1;
            }
            else if(slowcontrol::Variables::TaggerScalers->HasChanged())
            {
                tree.ExpLivetime()  = slowcontrol::Variables::Trigger->GetExpLivetime();
            }

            tree.cosThetaPi0COMS() = cos(getPi0COMS(tree.Tagg_E(), tree.photonSum()).Theta());
            tree.EMB_cosThetaPi0COMS() = cos(getPi0COMS(tree.EMB_Ebeam(), tree.EMB_photonSum()).Theta());
            tree.NCands() = data.Candidates.size();



            // write to tree:
            tree.Tree->Fill();

            if (FinalCuts())
                continue;
            FillStep("Final cuts");
            recSignal.Egamma()      = seenSignal.Egamma();
            recSignal.TaggerBin()   = seenSignal.TaggerBin();
            recSignal.CosThetaPi0() = seenSignal.CosThetaPi0();
            recSignal.Phi()         = seenSignal.Phi();
            recSignal.Theta()       = seenSignal.Theta();
            recSignal.Tree->Fill();
        }
    } // taggerHits
}

void singlePi0::Finish()
{
    if (!flag_mc)
        return;

    for ( long long en = 0 ; en < seenSignal.Tree->GetEntries() ; ++en )
    {
        seenSignal.Tree->GetEntry(en);
        hist_seen->Fill(seenSignal.TaggerBin(), seenSignal.CosThetaPi0());
    }

    for ( long long en = 0 ; en < recSignal.Tree->GetEntries() ; ++en )
    {
        recSignal.Tree->GetEntry(en);
        hist_rec->Fill(recSignal.TaggerBin(), recSignal.CosThetaPi0());
    }

    hist_efficiency->Add(hist_rec);
    hist_efficiency->Divide(hist_seen);
    //set to average for hadd later!!!
    hist_efficiency->SetBit(TH1D::kIsAverage);
}

void singlePi0::ShowResult()
{
    canvas("summary")
            << hist_steps
            << hist_channels
            << hist_channels_end
            << endc;

    canvas("eff")
            << drawoption("colz")
            << hist_efficiency
            << endc;
}

void singlePi0::PionProdTree::SetRaw(const utils::ProtonPhotonCombs::comb_t& selection)
{
    proton() = TSimpleParticle(*selection.Proton);
    photons() = TSimpleParticle::TransformParticleList(selection.Photons);

    photonSum() = selection.PhotonSum;
    IM2g()      = photonSum().M();
    IMproton_MM() = selection.MissingMass;
    DiscardedEk() = selection.DiscardedEk;

//    old, still needed?
//    proton_MM() = selection.MissingMass;
//    pg_copl    = selection.Copl_pg;
//    pMM_angle  = selection.Angle_pMM;
}



void singlePi0::PionProdTree::SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result)
{
    EMB_proton()     = TSimpleParticle(*(kF.GetFittedProton()));
    EMB_photons()    = TSimpleParticle::TransformParticleList(kF.GetFittedPhotons());
    EMB_photonSum()  = accumulate(EMB_photons().begin(),EMB_photons().end(),TLorentzVector(0,0,0,0));
    EMB_IM2g()       = EMB_photonSum().M();
    EMB_Ebeam()      = kF.GetFittedBeamE();
    EMB_iterations() = result.NIterations;
    EMB_prob()       = result.Probability;
    EMB_chi2()       = reducedChi2(result);
}



AUTO_REGISTER_PHYSICS(singlePi0)
