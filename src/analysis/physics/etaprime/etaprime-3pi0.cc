#include "etaprime-3pi0.h"
#include "utils/Combinatorics.h"
#include "base/std_ext/math.h"
#include "utils/ParticleTools.h"
#include "utils/ParticleID.h"
#include "base/ParticleType.h"
#include "base/ParticleTypeTree.h"

#include <algorithm>
#include <cassert>
#include <chrono>

#include "TTree.h"
#include "TCanvas.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;

Etap3pi0::Etap3pi0(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    signal_tree(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_3Pi0_6g)),
    reference_tree(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g)),
    bkg_tree(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)),
    fitterSig(utils::ParticleTools::GetProducedParticle(signal_tree),
              uncertModel, true,
              [] (ParticleTypeTree tree)
              {
                if(tree->Get() == ParticleTypeDatabase::EtaPrime)
                    return utils::TreeFitter::nodesetup_t{true};
                else
                    return utils::TreeFitter::nodesetup_t{};
               } ),
    fitterRef(utils::ParticleTools::GetProducedParticle(reference_tree),
              uncertModel, true,
              [] (ParticleTypeTree tree)
              {
                if(tree->Get() == ParticleTypeDatabase::EtaPrime)
                    return utils::TreeFitter::nodesetup_t{true};
                else
                    return utils::TreeFitter::nodesetup_t{};
              } ),
    kinFitterEMB(uncertModel,true)
{
        fitterSig.SetZVertexSigma(3);
        fitterRef.SetZVertexSigma(3);
        kinFitterEMB.SetZVertexSigma(3);

    promptrandom.AddPromptRange({-5,5});
    promptrandom.AddRandomRange({-20, -10});
    promptrandom.AddRandomRange({ 10,  20});

    tree = HistFac.makeTTree("tree");

    string cat("steps");
    AddHist1D(cat, "evcount",          "events after steps", "", "# events", BinSettings(5));

    cat = "channels";
    AddHist1D(cat,"mc_true",            "mc true for signal, ref, bkg", "", "#", BinSettings(3));
    AddHist1D(cat,"identified",         "IdentifiedChannels", "channel", "#", BinSettings(3));

    cat = "tagger";
    AddHist1D(cat,"tagHits",            "# Tagger hits", "# hits", "", BinSettings(15));


    vars.SetBranches(tree);

    auto SetupIntermediates = [] (const utils::TreeFitter& theFitter,vector<utils::TreeFitter::tree_t>& intermediates)
    {
        utils::TreeFitter::tree_t etaP = theFitter.GetTreeNode(ParticleTypeDatabase::EtaPrime);
        auto it = etaP->Daughters().begin();
        for (auto& pi: intermediates)
            pi = *it++;
    };

    SetupIntermediates(fitterSig,intermediatesTreeSig);
    SetupIntermediates(fitterRef,intermediatesTreeRef);
}

bool Etap3pi0::MakeMCProton(const TEventData& mcdata, TParticlePtr& proton)
{
    auto mctrue_particles = utils::ParticleTypeList::Make(mcdata.ParticleTree);
    const auto& protonlist = mctrue_particles.Get(ParticleTypeDatabase::Proton);
    if (protonlist.size() != 1)
        return false;
    proton = protonlist.at(0);
    return true;
}

double Etap3pi0::applyEnergyMomentumConservation(double EBeam, const ant::TParticleList& photons, const TParticlePtr& proton)
{
   APLCON::Result_t result;

   result = kinFitterEMB.DoFit(EBeam, proton, photons);
   vars.kinfitted.beamE = kinFitterEMB.GetFittedBeamE();
   vars.kinfitted.p = *kinFitterEMB.GetFittedProton();

   return result.ChiSquare;
}

void Etap3pi0::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);
    /// TODO:geo-cuts ??

    const auto& data   = event.Reconstructed();
    const auto& mcdata = event.MCTrue();

    TParticlePtr  proton;
    TParticleList photons;

    hists.at("steps").at("evcount")->Fill("1) totalEvts",1);

    // fill channels and set true particle ids if possible

    vars.truetype = -1;
    if ( mcdata.ParticleTree )
    {
        if (mcdata.ParticleTree->IsEqual(signal_tree, utils::ParticleTools::MatchByParticleName))
        {
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree).c_str(),1);
            vars.truetype = 0;
        }
        if (mcdata.ParticleTree->IsEqual(reference_tree, utils::ParticleTools::MatchByParticleName))
        {
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree).c_str(),1);
            vars.truetype = 1;
        }
        if (mcdata.ParticleTree->IsEqual(bkg_tree, utils::ParticleTools::MatchByParticleName))
        {
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree).c_str(),1);
            vars.truetype = 2;
        }
        vars.decayString = utils::ParticleTools::GetDecayString(mcdata.ParticleTree);
    }

    /// soon obsolete
    auto mctrue_particles = utils::ParticleTypeList::Make(mcdata.ParticleTree);
    const auto& mcprotons          = mctrue_particles.Get(ParticleTypeDatabase::Proton);

    if(!triggersimu.HasTriggered())
        return;
    hists.at("steps").at("evcount")->Fill("2) Triggered",1);
    vars.EsumCB = triggersimu.GetCBEnergySum();

    if ( data.Candidates.size() != 7)
        return;
    hists.at("steps").at("evcount")->Fill("3) 7 candidates",1);

    vars.protonTime = std_ext::NaN;
    for(const auto& cand : data.Candidates.get_iter())
    {
        if(cand->Detector & Detector_t::Type_t::TAPS)
        {
            if(!isfinite(vars.protonTime) || vars.protonTime < cand->Time)
            {
                vars.protonTime = cand->Time;
                proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand);
            }
        }
    }
    if (!proton)
        return;

    vars.proton = *proton;
    hists.at("steps").at("evcount")->Fill("4) proton in TAPS",1);
    vars.etaprimeCand = {};

    for ( const auto& cand: data.Candidates.get_iter() )
    {
       if ( cand.get_ptr() != proton->Candidate )
       {
         auto photon = make_shared<TParticle>(ParticleTypeDatabase::Photon, cand);
         vars.etaprimeCand += *photon;
         photons.emplace_back(move(photon));
       }
    }

    if (mcprotons.size() > 1)
        return;
    //debug:
    //hists.at("steps").at("evcount")->Fill("debug) <= 1 mc-true proton",1);
    if (mcprotons.size() == 1)
        vars.trueProton = *mcprotons.at(0);

    double CBAvgTime = triggersimu.GetRefTiming();
    if(!isfinite(CBAvgTime))
        return;
    //debug:
//    hists.at("steps").at("evcount")->Fill("debug) finite CBAvg-Time",1);

    // coplanarity
    vars.coplanarity = std_ext::radian_to_degree(vec2::Phi_mpi_pi(vars.proton.Phi() - vars.etaprimeCand.Phi() - M_PI ));
    if (fabs(vars.coplanarity) > phSettings.coplCut)
        return;
    hists.at("steps").at("evcount")->Fill((formatter() << "5) proton Coplanarity < " << phSettings.coplCut).str().c_str(),1);



    hists.at("tagger").at("tagHits")->Fill(data.TaggerHits.size());
    for(const TTaggerHit& t : data.TaggerHits )
    {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(t));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        vars.taggWeight    = promptrandom.FillWeight();
        vars.taggE         = t.PhotonEnergy;
        vars.taggCh        = unsigned(t.Channel);
        vars.taggTime      = t.Time;

        // cut on Tagger-Energy below eta' threshold
        if ( vars.taggE < phSettings.etaprimeThreshold )
            continue;
        hists.at("steps").at("evcount")->Fill("6) E_{#gamma} < E_{thresh}(#eta')",vars.taggWeight);

        LorentzVec beamPseudoParticle = t.GetPhotonBeam() + LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass());
        vars.MM = beamPseudoParticle - vars.etaprimeCand;

        assert( photons.size() == 6);

        // EMB - kinFit - cut
        vars.EMB_chi2 = applyEnergyMomentumConservation(t.PhotonEnergy,photons,proton);
        if ( vars.EMB_chi2 > phSettings.Chi2CutEMB )
            continue;
        hists.at("steps").at("evcount")->Fill((formatter() << "7) EMB-4C-KinFit: #chi^{2} < )" << phSettings.Chi2CutEMB).str().c_str(),
                                              vars.taggWeight);


        // IDF: ref & sig
//        MakeSignal(kinFitterEMB.GetFittedPhotons());
//        MakeReference(kinFitterEMB.GetFittedPhotons());
//        if (vars.chi2_sig < vars.chi2_ref)
//        {
//            if ( vars.chi2_sig < phSettings.Chi2CutSig )
//            {
//                vars.type = 0;
//                hists.at("steps").at("evcount")->Fill("8a) signal identified",vars.taggWeight);
//            }
//            else
//            {
//                vars.type = -1;
//                hists.at("steps").at("evcount")->Fill("8c) background identified",vars.taggWeight);

//            }
//        }
//        else
//        {
//            if (vars.chi2_ref < phSettings.Chi2CutRef )
//            {
//                vars.type = 1;
//                hists.at("steps").at("evcount")->Fill("8b) reference identified",vars.taggWeight);
//            }
//            else
//            {
//                vars.type = -1;
//                hists.at("steps").at("evcount")->Fill("8c) background identified",vars.taggWeight);
//            }
//        }
        tree->Fill();
    }
}

//void Etap3pi0::MakeSignal(const TParticleList& photonLeaves)
//{
//    fitterSig.SetPhotons(photonLeaves);
//    APLCON::Result_t result;
//    vars.chi2_sig = std::numeric_limits<double>::infinity();

//    while (fitterSig.NextFit(result))
//    {
//        if (result.Status != APLCON::Result_Status_t::Success )
//            continue;
//        if ( vars.chi2_sig < result.ChiSquare)
//            continue;
//        vars.chi2_sig      = result.ChiSquare;
//        vars.prob_sig      = result.Probability;
//        vars.iteration_sig = result.NIterations;

//        //signal is first, so fill 6g here...
//        vars.kinfitted.etaprimeCand = LorentzVec({0,0,0},0);
//        for (size_t i = 0 ; i < 3 ; ++i)
//        {
//            auto it_gamma =  intermediatesTreeSig[i]->Daughters().begin();
//            vars.kinfitted.gammasSig.at(2*i) = *((*it_gamma++)->Get().Leave->Particle);
//            vars.kinfitted.gammasSig.at((2*i)+1) = *((*it_gamma)->Get().Leave->Particle);
//            vars.kinfitted.intermediatesSig.at(i) =   vars.kinfitted.gammasSig.at(2*i)
//                                                    + vars.kinfitted.gammasSig.at((2*i)+1);
//            vars.kinfitted.etaprimeCand += vars.kinfitted.intermediatesSig.at(i);
//        }
//    }
//}

//void Etap3pi0::MakeReference(const TParticleList& photonLeaves)
//{
//    fitterRef.SetPhotons(photonLeaves);
//    APLCON::Result_t result;
//    vars.chi2_ref = std::numeric_limits<double>::infinity();

//    while (fitterRef.NextFit(result))
//    {
//        if (result.Status != APLCON::Result_Status_t::Success )
//            continue;
//        if ( vars.chi2_ref < result.ChiSquare)
//            continue;
//        vars.chi2_ref      = result.ChiSquare;
//        vars.prob_ref      = result.Probability;
//        vars.iteration_ref = result.NIterations;

//        for (size_t i = 0 ; i < 3 ; ++i)
//        {
//            auto it_gamma =  intermediatesTreeRef[i]->Daughters().begin();
//            vars.kinfitted.gammasRef.at(2*i) = *((*it_gamma++)->Get().Leave->Particle);
//            vars.kinfitted.gammasRef.at((2*i)+1) = *((*it_gamma)->Get().Leave->Particle);
//            vars.kinfitted.intermediatesRef.at(i) =   vars.kinfitted.gammasRef.at(2*i)
//                                                    + vars.kinfitted.gammasRef.at((2*i)+1);
//        }
//    }
//}

void Etap3pi0::Finish()
{
}

void Etap3pi0::ShowResult()
{
    vector<string> cats = { "steps"};
    for (auto& category: cats)
    {
        canvas c(category);
        for (auto h: hists.at(category))
            c << h.second;
        c << endc;
    }
}

void Etap3pi0::branches::SetBranches(TTree* tree)
{
    tree->Branch("proton",       &proton);
    tree->Branch("etaprimeCand", &etaprimeCand);
    tree->Branch("trueProton",   &trueProton);
    tree->Branch("MM",           &MM);

    tree->Branch("coplanarity", &coplanarity);

    tree->Branch("EsumCB", &EsumCB);

    tree->Branch("taggWeight", &taggWeight);
    tree->Branch("taggE",      &taggE);
    tree->Branch("taggCh",     &taggCh);
    tree->Branch("taggTime",   &taggTime);

    tree->Branch("EMB_chi2", &EMB_chi2);

    tree->Branch("chi2_ref",      &chi2_ref);
    tree->Branch("prob_ref",      &prob_ref);
    tree->Branch("iteration_ref", &iteration_ref);
    tree->Branch("status_ref",    &status_ref);

    tree->Branch("chi2_sig",      &chi2_sig);
    tree->Branch("prob_sig",      &prob_sig);
    tree->Branch("iteration_sig", &iteration_sig);
    tree->Branch("status_sig",    &status_sig);

    tree->Branch("type",     &type);
    tree->Branch("truetype", &truetype);

    tree->Branch("decayString", &decayString);

    tree->Branch("kf_beamE",      &kinfitted.beamE);
    tree->Branch("kf_inter_Sig",  &kinfitted.intermediatesSig);
    tree->Branch("kf_gammas_Sig", &kinfitted.gammasSig);
    tree->Branch("kf_gammas_Ref", &kinfitted.gammasRef);
    tree->Branch("kf_inter_Ref",  &kinfitted.intermediatesRef);
    tree->Branch("kf_6g",         &kinfitted.etaprimeCand);
    tree->Branch("kf_p",          &kinfitted.p);
}

void Etap3pi0::branches::FillKinfitBeamProton(double beamE, const ant::TParticlePtr& proton)
{
    kinfitted.beamE = beamE;
    kinfitted.p = *proton;
}

void Etap3pi0::AddHist1D(
        const std::string& category, const std::string& hname,
        const std::string& title,
        const std::string& xlabel, const std::string& ylabel,
        const BinSettings& bins
        )
{
    hists[category][hname] = HistFac.makeTH1D(title,xlabel,ylabel,bins,(category + string("_") + hname));
}

void Etap3pi0::AddHist2D(const string& category, const string& hname,
                         const string& title,
                         const string& xlabel, const string& ylabel,
                         const BinSettings& xbins, const BinSettings& ybins)
{
    hists[category][hname] = HistFac.makeTH2D(title,xlabel,ylabel,xbins,ybins,(category + string("_") + hname));
}

AUTO_REGISTER_PHYSICS(Etap3pi0)
