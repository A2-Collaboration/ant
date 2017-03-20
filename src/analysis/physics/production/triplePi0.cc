#include "triplePi0.h"

#include "expconfig/ExpConfig.h"

#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/vec/LorentzVec.h"
#include "base/Logger.h"

#include "utils/uncertainties/Interpolated.h"

#include "analysis/physics/Plotter.h"
#include "plot/CutTree.h"

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

    etaFitterBkg = fitterBkg.GetTreeNode(ParticleTypeDatabase::Eta);

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

    fitterBkg.SetIterationFilter([this] ()
    {
       const auto etaWindow = ParticleTypeDatabase::Eta.GetWindow(75);
       return  etaWindow.Contains(etaFitterBkg->Get().LVSum.M());
    });


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
                                           const tools::protonSelection_t& protonSelection)
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
    tree.CBAvgTime = triggersimu.GetRefTiming();

    auto bestProb_SIG  = 0.;
    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerHit));
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
            const auto selection =  tools::getProtonSelection(i_proton, data.Candidates,
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

                tree.SIG_photonVeto() = tools::getPhotonVetoEnergy(selection);
                tree.SIG_corrPhotonVeto() = tools::getPhotonVetoEnergy(selection,true);

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


        tree.ChargedClusterE() = tools::getChargedClusterE(data.Clusters);
        tree.ChargedCandidateE() = tools::getChargedCandidateE(data.Candidates);


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

void triplePi0::PionProdTree::SetRaw(const tools::protonSelection_t& selection)
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


using namespace ant::analysis::plot;

class triplePi0_PlotBase: public Plotter{

protected:
    TTree* t = nullptr;
    triplePi0::PionProdTree tree;
    // Plotter interface
public:
    triplePi0_PlotBase(const string& name, const WrapTFileInput& input, OptionsPtr opts):
        Plotter(name,input,opts)
    {
        if(!input.GetObject(triplePi0::treeAccessName(),t))
            throw Exception("Input TTree not found");

        if(!tree.Matches(t))
            throw runtime_error("Tree branches don't match");
        tree.LinkBranches(t);
    }

    virtual long long GetNumEntries() const override {return t->GetEntries();}
};

class triplePi0_Plot: public triplePi0_PlotBase {


protected:
    static const string data_name;
    static const double binScale;



    template<typename Hist_t>
    struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

        // Hist_t should have that type defined
        using Fill_t = typename Hist_t::Fill_t;



        MCTrue_Splitter(const HistogramFactory& histFac,
                        const cuttree::TreeInfo_t& treeInfo) :
            cuttree::StackedHists_t<Hist_t>(histFac, treeInfo)

        {
            using histstyle::Mod_t;

            // TODO: derive this from channel map
            this->GetHist(0, data_name, Mod_t::MakeDataPoints(kBlack));
            this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2));
            this->GetHist(2, "MainBkg",  Mod_t::MakeLine(kGreen, 2));
            // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
            // see also Fill()
            this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
            this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));

        }

        void Fill(const Fill_t& f) {

            const unsigned mctrue = f.Tree.MCTrue;

            using histstyle::Mod_t;

            auto get_bkg_name = [] (const unsigned mctrue) {
                return physics::triplePi0::getOtherChannelNames(mctrue); //(int(mctrue));
            };

            using histstyle::Mod_t;

            const Hist_t& hist = mctrue<9 ? this->GetHist(mctrue) :
                                            this->GetHist(mctrue,
                                                           get_bkg_name(mctrue),
                                                           Mod_t::MakeLine(histstyle::color_t::GetLight(mctrue-10), 1, kGray+1)
                                                           );


            hist.Fill(f);

            // handle MC_all and MC_bkg
            if(mctrue>0) {
                this->GetHist(3).Fill(f);
                if(mctrue >= 9 || mctrue == 2)
                    this->GetHist(4).Fill(f);
            }
        }
    };

    struct TriplePi0Hist_t {

        using Tree_t = physics::triplePi0::PionProdTree;

        struct Fill_t {
            const Tree_t& Tree;

            Fill_t(const Tree_t& t) : Tree(t) {}

            double TaggW() const {
                return Tree.Tagg_W;
            }

            vector<TLorentzVector> get2G(const vector<TLorentzVector>& photons) const
            {
                vector<TLorentzVector> acc;
                const auto& permutation(Tree.SIG_combination());
                for (size_t i = 0 ; i < permutation.size(); i+=2)
                {
                    TLorentzVector gg(photons.at(permutation.at(i)));
                    gg += photons.at(permutation.at(i+1));
                    acc.emplace_back(gg);
                }
                return acc;
            }
        };

        template <typename Hist>
        using fillfunc_t = std::function<void(Hist*, const Fill_t&)>;

        template <typename Hist>
        struct HistFiller_t {
            fillfunc_t<Hist> func;
            Hist* h;
            HistFiller_t(Hist* hist, fillfunc_t<Hist> f): func(f), h(hist) {}
            void Fill(const Fill_t& data) const {
                func(h, data);
            }
        };

        template <typename Hist>
        struct HistMgr : std::list<HistFiller_t<Hist>> {

            using list<HistFiller_t<Hist>>::list;

            void Fill(const Fill_t& data) const {
                for(auto& h : *this) {
                    h.Fill(data);
                }
            }
        };

        static BinSettings Bins(const unsigned bins, const double min, const double max) {
            return BinSettings(unsigned(bins*binScale), min, max);
        }

        HistMgr<TH1D> h1;
        HistMgr<TH2D> h2;

        const BinSettings Ebins    = Bins(1000, 0, 1000);

        const BinSettings Chi2Bins = Bins(250, 0,   25);
        const BinSettings probbins = Bins(250, 0,   1);

        const BinSettings IMbins       = Bins(1000,  200, 1100);
        const BinSettings IMProtonBins = Bins(1000,  600, 1200);
        const BinSettings IM2g         = Bins(1000,    0,  360);

        const BinSettings pThetaBins = Bins( 200,  0,   80);
        const BinSettings pEbins     = Bins( 350,  0, 1200);

        HistogramFactory HistFac;

        void AddTH1(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name, fillfunc_t<TH1D> f) {
            h1.emplace_back(HistFiller_t<TH1D>(
                                HistFac.makeTH1D(title, xlabel, ylabel, bins, name),f));
        }

        void AddTH2(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings& ybins, const string &name, fillfunc_t<TH2D> f) {
            h2.emplace_back(HistFiller_t<TH2D>(
                                HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name),f));
        }

        TriplePi0Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf)
        {
            AddTH1("TreeFit Probability",      "probability",             "",       probbins,   "TreeFitProb",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.SIG_prob, f.TaggW());
            });

            AddTH1("6#gamma IM","6#gamma IM [MeV]", "", IMbins,"IM_6g",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IM6g, f.TaggW());
            });

    //        AddTH1("Sig && Bkg", "6#gammaa IM [MeV]", "",IMbins,"IM_6g_correct",
    //               [] (TH1D* h, const Fill_t& f)
    //        {
    //            auto correctF = f.TaggW();
    //            if (!(f.Tree.SIG_combination().size() == 0 || f.Tree.BKG_combination().size() == 0 ))
    //                for ( auto i = 0u ; i < f.Tree.SIG_combination().size() ; ++i)
    //                    if (f.Tree.SIG_combination().at(i) != f.Tree.BKG_combination().at(i))
    //                    {
    //                        correctF = 0.0;
    //                        break;
    //                    }
    //            h->Fill(f.Tree.EMB_IM6g, correctF);
    //        });

            AddTH1("6#gamma IM fitted","6#gamma IM [MeV]", "", IMbins,"IM_6g_fit",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_IM6g, f.TaggW());
            });

            AddTH1("tree fitted 3#pi^{0}","IM_{3#pi^{0}} [MeV]","",IMbins,"3pi0im",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.SIG_IM3Pi0,f.TaggW());
            });

            AddTH1("2g MM SIG combination","MM_{2#gamma} [MeV]","",IM2g,"combSig2g",
                   [] (TH1D* h, const Fill_t& f)
            {
                const auto gammas = f.get2G(f.Tree.EMB_photons());
                for( const auto& m: gammas)
                    h->Fill(m.M(),f.TaggW());
            });


            AddTH1("MM proton","MM_{proton} [MeV]", "", IMProtonBins, "IM_p",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.proton_MM().M(), f.TaggW());
            });

            AddTH1("MM proton fitted","MM_{proton} [MeV]", "", IMProtonBins, "IM_p_fit",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton_MM().M(),f.TaggW());
            });

            AddTH1("MM pions", "IM_{2#gamma} [MeV]","", IM2g,"IM_pions",
                   [] (TH1D* h, const Fill_t& f)
            {
                for ( const auto& pion: f.Tree.SIG_pions())
                    h->Fill(pion.M(),f.TaggW());
            });

            AddTH1("Proton_MM_Angle", "Angle [#circ]","", Bins(200,0,40),"MM_pAngle",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.pMM_angle,f.TaggW());
            });

            AddTH1("CB_ESum", "EsumCB [MeV]","", Bins(300,500,1900),"CBESUM",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.CBESum, f.TaggW());
            });

            AddTH2("Fitted Proton","E^{kin}_{p} [MeV]","#theta_{p} [#circ]",pEbins,pThetaBins,"pThetaVsE",
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton().E() - ParticleTypeDatabase::Proton.Mass(), std_ext::radian_to_degree(f.Tree.EMB_proton().Theta()), f.TaggW());
            });

            AddTH2("Resonance Search 1","m(p #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",Bins(300,  900, 1900),Bins(300,    0, 1000),"ppi0_2pi0",
                   [] (TH2D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                for (auto i = 0u; i < pions.size() ; ++i)
                {
                    const auto N    = pions.at(i) + proton;
                    LorentzVec pipi({0,0,0},0);
                    for (auto j = 0u; j < pions.size() ; ++j)
                        if ( j != i )
                            pipi += pions.at(j);

                    h->Fill(N.M(),pipi.M(),f.TaggW());
                }
            });

            AddTH1("Resonance Search 1","m(p #pi^{0}) [MeV]","",Bins(300,  900, 1900),"ppi0",
                   [] (TH1D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                if (pions.size() == 3)
                {
                    for(auto i = 0u ; i < 3 ; ++i)
                    {
                        const auto N    = pions.at(i) + f.Tree.EMB_proton();
                        h->Fill(N.M(),f.TaggW());
                    }
                } else { LOG(INFO) << pions.size() ;}
            });

            AddTH2("Resonance Search 2","m(2 #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",Bins(300,  0, 1000),Bins(300,    0, 1000),"2pi0_2pi0",
                   [] (TH2D* h, const Fill_t& f)
            {
                const vector<pair<size_t,size_t>> combinations = { { 0 , 1 } , { 0 , 2 } , { 1 , 2 } };
                const auto pions = f.Tree.SIG_pions();
                if (pions.size() == 0)
                    return;

                for ( size_t i = 0 ; i < 3 ; ++i)
                    for ( size_t j = 0 ; j < 3 ; ++j)
                    {
                        if ( i == j )
                            continue;
                        const auto ppM2  =(pions.at(combinations.at(i).first) + pions.at(combinations.at(i).second)).M();
                        const auto ppM1  =(pions.at(combinations.at(j).first) + pions.at(combinations.at(j).second)).M();
                        h->Fill(ppM2,ppM1,f.TaggW());
                    }
            });

            AddTH1("Resonance Search 3","m(p 2 #pi^{0}) [MeV]","",Bins(1000, 1000, 2000),"p2pi0",
                   [] (TH1* h, const Fill_t& f)
            {
    //            const auto pions = f.get2G(f.Tree.EMB_photons());
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();
                for(auto comb=utils::makeCombination(pions, 2); !comb.done(); ++comb)
                {
                    const auto N  = comb.at(0) + comb.at(1) + f.Tree.EMB_proton();
                    h->Fill(N.M(),f.TaggW());
                }
            });

        }

        void Fill(const Fill_t& f) const {
            h1.Fill(f);
            h2.Fill(f);
        }

        std::vector<TH1*> GetHists() const {
            vector<TH1*> v;
            v.reserve(h1.size()+h2.size());
            for(auto& e : h1) {
                v.emplace_back(e.h);
            }
            for(auto& e: h2) {
                v.emplace_back(e.h);
            }
            return v;
        }


//        static TCutG* makeDalitzCut() {
//            TCutG* c = new TCutG("DalitzCut", 3);
//            c->SetPoint(0, 0.0,  0.2);
//            c->SetPoint(1, -.22, -.11);
//            c->SetPoint(2,  .22, -.11);
//            return c;
//        }

//        static TCutG* dalitzCut;

        struct TreeCuts {

            static bool KinFitProb(const Fill_t& f) noexcept {
                return     f.Tree.EMB_prob >  0.1;
            }
            static bool proton_MM(const Fill_t& f) noexcept {
                const auto width = 180.0;
                const auto mmpm = f.Tree.proton_MM().M();
                return (938.3 - width < mmpm && mmpm < 938.3 + width);
            }
        };

        // Sig and Ref channel share some cuts...
        static cuttree::Cuts_t<Fill_t> GetCuts() {

            using cuttree::MultiCut_t;

            cuttree::Cuts_t<Fill_t> cuts;

            const cuttree::Cut_t<Fill_t> ignore({"ignore", [](const Fill_t&){ return true; }});


            cuts.emplace_back(MultiCut_t<Fill_t>{
                                 { "EMB_prob > 0.1", [](const Fill_t& f)
                                   {
                                       return TreeCuts::proton_MM(f);
                                   }
                                 },
                                  ignore
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"SIG_prob > 0.1", [](const Fill_t& f)
                                   {
                                       return f.Tree.SIG_prob > 0.1;
                                   }
                                  }
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"BKG_prob < 0.1", [](const Fill_t& f)
                                   {
                                       return f.Tree.BKG_prob < 0.1;
                                   }
                                  },
                                  {
                                      "IM 6g >  600 MeV", [](const Fill_t& f)
                                      {
                                          return f.Tree.SIG_IM3Pi0 > 600;
                                      }
                                  }
                              });



            return cuts;
        }

    };

    plot::cuttree::Tree_t<MCTrue_Splitter<TriplePi0Hist_t>> signal_hists;

    // Plotter interface
public:

    triplePi0_Plot(const string& name, const WrapTFileInput& input, OptionsPtr opts):
        triplePi0_PlotBase(name,input,opts)
    {
        signal_hists = cuttree::Make<MCTrue_Splitter<TriplePi0Hist_t>>(HistFac);
    }


    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);
        cuttree::Fill<MCTrue_Splitter<TriplePi0Hist_t>>(signal_hists, {tree});
    }

    virtual void Finish() override{}
    virtual void ShowResult() override{}

    virtual ~triplePi0_Plot(){}

};

class triplePi0_Test: public triplePi0_PlotBase{

protected:

    TH1D* m3pi0  = nullptr;

    unsigned nBins  = 200u;



    bool testCuts() const
    {
        return (
                    true                      &&
                    tree.SIG_prob   < 0.1     &&
                    tree.SIG_IM3Pi0 < 600.0
                );
    }

    // Plotter interface
public:
    triplePi0_Test(const string& name, const WrapTFileInput& input, OptionsPtr opts):
        triplePi0_PlotBase (name,input,opts)
    {
        m3pi0  = HistFac.makeTH1D("m(3#pi^{0})","m(3#pi^0) [MeV]","#",      BinSettings(nBins,400,1100));
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



const string triplePi0_Plot::data_name = "Data";
const double triplePi0_Plot::binScale  = 1.0;


AUTO_REGISTER_PHYSICS(triplePi0)
AUTO_REGISTER_PLOTTER(triplePi0_Plot)
AUTO_REGISTER_PLOTTER(triplePi0_Test)

