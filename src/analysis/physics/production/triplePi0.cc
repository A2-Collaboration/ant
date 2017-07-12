#include "triplePi0.h"

#include "expconfig/ExpConfig.h"

#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/vec/LorentzVec.h"
#include "base/Logger.h"

#include "utils/uncertainties/Interpolated.h"

#include "analysis/physics/Plotter.h"
#include "plot/CutTree.h"

#include "slowcontrol/SlowControlVariables.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;



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



auto getEgamma = [] (const TParticleTree_t& tree)
{
    return tree->Get()->Ek();
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


//auto getLorentzSumUnfitted = [](const vector<utils::TreeFitter::tree_t>& nodes)
//{
//    vector<TLorentzVector> acc;
//    for ( const auto& node: nodes)
//    {
//        LorentzVec temp({0,0,0},0);
//        for ( const auto& ph: node->Daughters())
//        {
//            temp+=(*(ph->Get().Leave->Particle));
//        }
//        acc.push_back(temp);
//    }
//    return acc;
//};

//auto getTLorentz = [] (const utils::TreeFitter::tree_t& node)
//{
//    return node->Get().LVSum;
//};

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
    fitterEMB(                              uncertModelData, true ),
    fitterSig(signal.DecayTree,                uncertModelData, true )
//    fitterSigmaPlus(sigmaBackground.DecayTree, uncertModel, true )
{

    fitterSig.SetZVertexSigma(phSettings.fitter_ZVertex);
//    fitterSigmaPlus.SetZVertexSigma(phSettings.fitter_ZVertex);
    fitterEMB.SetZVertexSigma(phSettings.fitter_ZVertex);


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


//    pionsFitterSigmaPlus = fitterSigmaPlus.GetTreeNodes(ParticleTypeDatabase::Pi0);

//    kaonFitterSigmaPlus  = fitterSigmaPlus.GetTreeNode(ParticleTypeDatabase::K0s);
//    sigmaFitterSigmaPlus = fitterSigmaPlus.GetTreeNode(ParticleTypeDatabase::SigmaPlus);

//    // be lazy and catch complete class...
//    fitterSigmaPlus.SetIterationFilter([this] () {
//        const auto sigmaPlus_cut = ParticleTypeDatabase::SigmaPlus.GetWindow(200);
//        const auto K0s_cut = ParticleTypeDatabase::K0s.GetWindow(100);
//        auto ok = sigmaPlus_cut.Contains(sigmaFitterSigmaPlus->Get().LVSum.M()) &&
//                  K0s_cut.Contains(kaonFitterSigmaPlus->Get().LVSum.M());
//        return ok;
//    });



    promptrandom.AddPromptRange(phSettings.Range_Prompt);
    for ( const auto& range: phSettings.Ranges_Random)
        promptrandom.AddRandomRange(range);

    hist_steps          = HistFac.makeTH1D("steps","","# evts.",BinSettings(1,0,0),"steps");
    hist_channels       = HistFac.makeTH1D("channels","","# evts.",BinSettings(1,0,0),"channels");
    hist_channels_end   = HistFac.makeTH1D("channel-selected","","# evts.",BinSettings(1,0,0),"channels_end");

    hist_neutrals_channels
            = HistFac.makeTH2D("# neutral candidates","","# neutrals",BinSettings(1,0,0),BinSettings(15),"channels_neutrals");


    tree.CreateBranches(HistFac.makeTTree(phSettings.Tree_Name));
    seenSignal.CreateBranches(HistFac.makeTTree(seenSignal.treeName()));
    recSignal.CreateBranches(HistFac.makeTTree(recSignal.treeName()));

    if (!flag_mc)
    {
        slowcontrol::Variables::TaggerScalers->Request();
        slowcontrol::Variables::Trigger->Request();
        slowcontrol::Variables::TaggEff->Request();
    }
    fitterEMB.SetUncertaintyModel(flag_mc ? uncertModelMC : uncertModelData);
    fitterSig.SetUncertaintyModel(flag_mc ? uncertModelMC : uncertModelData);

}

triplePi0::fitRatings_t applyTreeFit(utils::TreeFitter& fitter,
                                     const std::vector<utils::TreeFitter::tree_t>& intermediates,
                                     const utils::ProtonPhotonCombs::comb_t& protonSelection,
                                     const double Ebeam)
{

    fitter.PrepareFits(Ebeam,
                       protonSelection.Proton,
                       protonSelection.Photons);
    APLCON::Result_t result;
    auto best_prob = std_ext::NaN;
    triplePi0::fitRatings_t fr(0,0,0,false,
                               {},
                               {},{});
    while(fitter.NextFit(result))
        if (   (result.Status    == APLCON::Result_Status_t::Success)
               && (std_ext::copy_if_greater(best_prob,result.Probability)))
        {

            fr = triplePi0::fitRatings_t(best_prob,reducedChi2(result),result.NIterations,
                                         result.Status == APLCON::Result_Status_t::Success,
                                         TSimpleParticle(*fitter.GetFittedProton()),
                                         getLorentzSumFitted(intermediates),
                                         getTreeFitPhotonIndices(protonSelection.Photons,fitter));
        }

    return fr;
}

void triplePi0::ProcessEvent(const ant::TEvent& event, manager_t&)
{
    const auto& data   = event.Reconstructed();

    triggersimu.ProcessEvent(event);


    FillStep("seen");

    tree.CBESum = triggersimu.GetCBEnergySum();
    // check if mc-flag ist set properly:
    if ( flag_mc != data.ID.isSet(TID::Flags_t::MC))
        throw runtime_error("provided mc flag does not match input files!");



//    const auto& mcTrue       = event.MCTrue();
    const auto& particleTree = event.MCTrue().ParticleTree;
    //===================== TreeMatching   ====================================================
    string trueChannel = "data";
    tree.MCTrue = phSettings.Index_Data;
    if (flag_mc)
    {
        tree.MCTrue() = phSettings.Index_brokenTree;
        trueChannel = "no Pluto tree";
    }
    if (particleTree)
    {

        if (particleTree->IsEqual(signal.DecayTree,utils::ParticleTools::MatchByParticleName))
        {
            const auto& taggerhits = event.MCTrue().TaggerHits;

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
                seenSignal.Tree->Fill();
            }
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
                index++;
                if (particleTree->IsEqual(otherChannel.DecayTree,utils::ParticleTools::MatchByParticleName))
                {
                    tree.MCTrue = index;
                    trueChannel = otherChannel.Name;
                    found = true;
                }
            }
            if (!found)
            {
                tree.MCTrue = phSettings.Index_Offset;
                trueChannel = utils::ParticleTools::GetDecayString(particleTree) + ": unknown";
            }
        }

    }
    hist_channels->Fill(trueChannel.c_str(),1);

    //simulate cb-esum-trigger
    if (!triggersimu.HasTriggered())
        return;
    FillStep("Triggered");

    if (tools::cutOn("N_{cands}",phSettings.Cut_NCands,data.Candidates.size(),hist_steps)) return;


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

        tree.Tagg_Ch  = static_cast<unsigned>(taggerHit.Channel);
        tree.Tagg_E   = taggerHit.PhotonEnergy;
        tree.Tagg_W   = promptrandom.FillWeight();

        {
            const auto taggEff = slowcontrol::Variables::TaggEff->Get(taggerHit.Channel);
            tree.Tagg_Eff      = taggEff.Value;
            tree.Tagg_EffErr  = taggEff.Error;
        }


        auto selections =  proton_photons()
                           .FilterMult(phSettings.nPhotons,100)
                           .FilterMM(taggerHit, phSettings.Cut_MM);
        if (selections.empty())
        {
            FillStep("No combs left");
            continue;
        }

        auto bestFitProb = 0.0;
        auto bestFound   = false;
        for ( const auto& selection: selections)
        {
            /// constraint fits
            const auto EMB_result = fitterEMB.DoFit(taggerHit.PhotonEnergy, selection.Proton, selection.Photons);
            if (!(EMB_result.Status == APLCON::Result_Status_t::Success))
                continue;
            FillStep("EMB-fit success");
            if (tools::cutOn("EMB-prob",phSettings.Cut_EMB_prob,EMB_result.Probability,hist_steps)) continue;
            const auto sigFitRatings = applyTreeFit(fitterSig,pionsFitterSig,selection,taggerHit.PhotonEnergy);
            if (!(sigFitRatings.FitOk))
                continue;
            FillStep("Tree-Fit succesful");


            ///status:
            const auto prob = phSettings.selType == settings_t::selectOn::kinFit ? EMB_result.Probability
                                                          : sigFitRatings.Prob;


            if ( prob > bestFitProb )
            {
                bestFound = true;
                bestFitProb = prob;

                tree.SetRaw(selection);
                tree.SetEMB(fitterEMB,EMB_result);
                tree.SetSIG(sigFitRatings);



            } // endif best SIG - treefit

        } // proton - candidate - loop


        if (bestFound)
        {
            FillStep("p identified");
            tree.ChargedClusterE() = tools::getChargedClusterE(data.Clusters);
            tree.ChargedCandidateE() = tools::getCandidateVetoE(data.Candidates);

            const auto neutralCands = tools::getNeutral(data,phSettings.vetoThreshE);

            tree.Neutrals() = neutralCands.size();

            if (flag_mc)
            {
                tree.ExpLivetime() = 1;
            }
            else if(slowcontrol::Variables::TaggerScalers->HasChanged())
            {
                tree.ExpLivetime()  = slowcontrol::Variables::Trigger->GetExpLivetime();
            }
            tree.NCands() = data.Candidates.size();

            recSignal.Egamma()      = seenSignal.Egamma();
            recSignal.TaggerBin()   = seenSignal.TaggerBin();

            hist_channels_end->Fill(trueChannel.c_str(),1);
            hist_neutrals_channels->Fill(trueChannel.c_str(),neutralCands.size(),1);
            recSignal.Tree->Fill();
            tree.Tree->Fill();
        }

    } // taggerHits - loop
}

void triplePi0::ShowResult()
{
    const auto colz = drawoption("colz");


    canvas("summary")
            << hist_steps
            << hist_channels
            << hist_channels_end
            << TTree_drawable(tree.Tree,"IM6g")
            << endc;

    canvas("channels")
            << colz
            << hist_neutrals_channels
            << endc;
}

void triplePi0::PionProdTree::SetRaw(const utils::ProtonPhotonCombs::comb_t& selection)
{
    proton()      = TSimpleParticle(*selection.Proton);
    photons()     = TSimpleParticle::TransformParticleList(selection.Photons);

    ProtonVetoE()  = selection.Proton->Candidate->VetoEnergy;
    PionPIDVetoE() = 0.;
    for (const auto& g: selection.Photons)
    {
        for (const auto& c: g->Candidate->Clusters)
        {
            if (c.DetectorType == Detector_t::Type_t::PID)
                PionPIDVetoE += c.Energy;
        }
    }


    photonSum()   = selection.PhotonSum;
    IM6g()        = photonSum().M();
    proton_MM()   = selection.MissingMass;
    DiscardedEk() = selection.DiscardedEk;
}



void triplePi0::PionProdTree::SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result)
{

    const auto fittedPhotons = kF.GetFittedPhotons();
    const auto phE           = kF.GetFittedBeamE();

    EMB_proton     = TSimpleParticle(*(kF.GetFittedProton()));
    EMB_photons    = TSimpleParticle::TransformParticleList(fittedPhotons);
    EMB_photonSum  = accumulate(EMB_photons().begin(),EMB_photons().end(),LorentzVec({0,0,0},0));
    EMB_IM6g       = EMB_photonSum().M();
    EMB_Ebeam      = phE;
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
    SIG_IM6g        = accumulate(fitRating.Intermediates.begin(),
                                 fitRating.Intermediates.end(),
                                 LorentzVec({0,0,0},0)).M();
    SIG_proton      = fitRating.Proton;
    SIG_combination = fitRating.PhotonCombination;
}


using namespace ant::analysis::plot;

using triplePi0_PlotBase = TreePlotterBase_t<triplePi0::PionProdTree>;



class triplePi0_Test: public triplePi0_PlotBase{

protected:

    TH1D* m3pi0  = nullptr;

    unsigned nBins  = 200u;



    bool testCuts() const
    {
        return (
                    true                      &&
                    tree.SIG_prob() < 0.1     &&
                    tree.SIG_IM6g() < 600.0
                );
    }

    // Plotter interface
public:
    triplePi0_Test(const string& name, const WrapTFileInput& input, OptionsPtr opts):
        triplePi0_PlotBase (name,input,opts)
    {
        m3pi0  = HistFac.makeTH1D("m(3#pi^{0})","m(3#pi^0) [MeV]","#",
                                  BinSettings(nBins,400,1100));
    }

    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);

        if (testCuts()) return;


        m3pi0->Fill(tree.IM6g);
    }
    virtual void ShowResult() override
    {
        canvas("view")
                << m3pi0
                << endc;
    }
};



AUTO_REGISTER_PHYSICS(triplePi0)
AUTO_REGISTER_PLOTTER(triplePi0_Test)

// code snippet for SIGMA treefit:
//                tree.SetBKG(applyTreeFit(fitterBkg,pionsFitterBkg,selection));
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
