#include "etaprime-3pi0.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "base/std_ext/math.h"
#include "utils/particle_tools.h"
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
    fitterSig("fitterSig",utils::ParticleTools::GetProducedParticle(signal_tree),
              [] (ParticleTypeTree tree)
              {
                if(tree->Get() == ParticleTypeDatabase::EtaPrime)
                    return utils::TreeFitter::nodesetup_t{1.0, true};
                else
                    return utils::TreeFitter::nodesetup_t{};
               } ),
    fitterRef("fitterRef",utils::ParticleTools::GetProducedParticle(reference_tree),
              [] (ParticleTypeTree tree)
              {
                if(tree->Get() == ParticleTypeDatabase::EtaPrime)
                    return utils::TreeFitter::nodesetup_t{1.0, true};
                else
                    return utils::TreeFitter::nodesetup_t{};
              } ),
    kinFitterEMB(GetName(), 6)
{
    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

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


    kinFitterEMB.SetupBranches(tree, "kinFitEMB");
    kinFitterEMB.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
    fitterSig.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
    fitterRef.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    vars.SetBranches(tree);


    utils::TreeFitter::tree_t etaP = fitterSig.GetTreeNode(ParticleTypeDatabase::EtaPrime);
    auto it = etaP->Daughters().begin();
    pi1Sig = *it++;
    pi2Sig = *it++;
    pi3Sig = *it;

}

bool Etap3pi0::MakeMCProton(const TEventData& mcdata, TParticlePtr& proton)
{
   const auto& protonlist = mcdata.Particles.Get(ParticleTypeDatabase::Proton);
   if (protonlist.size() != 1)
       return false;
   proton = protonlist.at(0);
   return true;
}

double Etap3pi0::getEnergyMomentumConservation(double EBeam, const ant::TParticleList& photons, const TParticlePtr& proton)
{
   APLCON::Result_t result;

   kinFitterEMB.SetEgammaBeam(EBeam);
   kinFitterEMB.SetPhotons(photons);
   kinFitterEMB.SetProton(proton);

   result = kinFitterEMB.DoFit();

   vars.FillKinFit( kinFitterEMB.GetFittedBeamE(),
                    kinFitterEMB.GetFittedPhotons(),
                    kinFitterEMB.GetFittedProton());

   return result.ChiSquare;
}

void Etap3pi0::ProcessEvent(const TEvent& event, manager_t&)
{
    /// TODO:geo-cuts ??

    const auto& data   = *event.Reconstructed;
    const auto& mcdata = *event.MCTrue;

    // maybe this????
    //TParticleList intermediate_SIG;
    //TParticleList intermediate_REF;

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
    const auto& mcprotons          = mcdata.Particles.Get(ParticleTypeDatabase::Proton);

    if(data.Trigger.CBEnergySum < phSettings.EsumCB)
        return;
    hists.at("steps").at("evcount")->Fill("2) CB-Energy-Sum",1);
    vars.EsumCB = data.Trigger.CBEnergySum;

    if ( data.Candidates.size() != 7)
        return;
    hists.at("steps").at("evcount")->Fill("3) 7 cands",1);

    vars.protonTime = std_ext::NaN;
    for(const auto& cand : data.Candidates)
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

    for ( const auto& cand: data.Candidates )
    {
       if ( cand != proton->Candidate )
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

    double CBAvgTime = event.Reconstructed->Trigger.CBTiming;
    if(!isfinite(CBAvgTime))
        return;
    //debug:
//    hists.at("steps").at("evcount")->Fill("debug) finite CBAvg-Time",1);

    // coplanarity
    vars.coplanarity = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(vars.proton.Phi() - vars.etaprimeCand.Phi() - M_PI ));
    if (fabs(vars.coplanarity) > phSettings.coplCut)
        return;
    hists.at("steps").at("evcount")->Fill("5) proton Coplanarity",1);



    hists.at("tagger").at("tagHits")->Fill(data.TaggerHits.size());
    for(const TTaggerHit& t : data.TaggerHits )
    {
        promptrandom.SetTaggerHit(t.Time - CBAvgTime);
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

        TLorentzVector beamPseudoParticle = t.GetPhotonBeam() + TLorentzVector(0,0,0, ParticleTypeDatabase::Proton.Mass());
        vars.MM = beamPseudoParticle - vars.etaprimeCand;

        assert( photons.size() == 6);

        // EMB - kinFit - cut
        vars.EMB_chi2 = getEnergyMomentumConservation(t.PhotonEnergy,photons,proton);
        if ( vars.EMB_chi2 > phSettings.fourConstrainChi2Cut )
            continue;
        hists.at("steps").at("evcount")->Fill("7) EMB-4C-KinFit: #chi^{2} < 40)",vars.taggWeight);


        // IDF: ref & sig
        MakeSignal(photons);
        MakeReference(photons);
        if (vars.chi2_sig < vars.chi2_ref)
        {
            if ( vars.chi2_sig < phSettings.fourConstrainChi2Cut )
            {
                vars.type = 0;
                hists.at("steps").at("evcount")->Fill("8a) signal identified",vars.taggWeight);
            }
            else
            {
                vars.type = -1;
                hists.at("steps").at("evcount")->Fill("8c) background identified",vars.taggWeight);

            }
        }
        else
        {
            if (vars.chi2_ref < phSettings.fourConstrainChi2Cut )
            {
                vars.type = 1;
                hists.at("steps").at("evcount")->Fill("8b) reference identified",vars.taggWeight);
            }
            else
            {
                vars.type = -1;
                hists.at("steps").at("evcount")->Fill("8c) background identified",vars.taggWeight);
            }
        }
        tree->Fill();
    }
}

void Etap3pi0::MakeSignal(const TParticleList& photonLeaves)
{
    fitterSig.SetLeaves(photonLeaves);
    APLCON::Result_t result;
    vars.chi2_sig = std::numeric_limits<double>::infinity();

    while (fitterSig.NextFit(result))
    {
        if (result.Status != APLCON::Result_Status_t::Success )
            continue;
        if ( vars.chi2_sig < result.ChiSquare)
            continue;
        vars.chi2_sig      = result.ChiSquare;
        vars.prob_sig      = result.Probability;
        vars.iteration_sig = result.NIterations;
        auto it_gamma =  pi1Sig->Daughters().begin();
        vars.kinfittedSig.gammas.at(0) = *((*it_gamma++)->Get().Leave->Particle);
        vars.kinfittedSig.gammas.at(1) = *((*it_gamma)->Get().Leave->Particle);
        vars.kinfittedSig.intermediates.at(0) = vars.kinfittedSig.gammas.at(0) + vars.kinfittedSig.gammas.at(1);
        it_gamma =  pi2Sig->Daughters().begin();
        vars.kinfittedSig.gammas.at(2) = *((*it_gamma++)->Get().Leave->Particle);
        vars.kinfittedSig.gammas.at(3) = *((*it_gamma)->Get().Leave->Particle);
        vars.kinfittedSig.intermediates.at(1) = vars.kinfittedSig.gammas.at(2) + vars.kinfittedSig.gammas.at(3);
        it_gamma =  pi3Sig->Daughters().begin();
        vars.kinfittedSig.gammas.at(4) = *((*it_gamma++)->Get().Leave->Particle);
        vars.kinfittedSig.gammas.at(5) = *((*it_gamma)->Get().Leave->Particle);
        vars.kinfittedSig.intermediates.at(2) = vars.kinfittedSig.gammas.at(4) + vars.kinfittedSig.gammas.at(5);
    }
}

void Etap3pi0::MakeReference(const TParticleList& photonLeaves)
{
    fitterRef.SetLeaves(photonLeaves);
    APLCON::Result_t result;
    vars.chi2_ref = std::numeric_limits<double>::infinity();

    while (fitterRef.NextFit(result))
    {
        if (result.Status != APLCON::Result_Status_t::Success )
            continue;
        if ( vars.chi2_ref < result.ChiSquare)
            continue;
        vars.chi2_ref      = result.ChiSquare;
        vars.prob_ref      = result.Probability;
        vars.iteration_ref = result.NIterations;
    }
}

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
    cout << "Setting up branches" << endl;
    tree->Branch("proton", &proton);

    tree->Branch("etaprimeCand", &etaprimeCand);

    tree->Branch("fittedProton",&fittedProton);

    tree->Branch("trueProton", &trueProton);

    tree->Branch("MM", &MM);

    tree->Branch("coplanarity", &coplanarity);
    tree->Branch("EsumCB", &EsumCB);

    tree->Branch("taggWeight", &taggWeight);
    tree->Branch("taggE", &taggE);
    tree->Branch("taggCh", &taggCh);
    tree->Branch("taggTime", &taggTime);

    tree->Branch("EMB_chi2",&EMB_chi2);

    tree->Branch("pi0s", &pi0);
    tree->Branch("pi0_chi2[3]", pi0_chi2, "pi0_chi2[3]/D");
    tree->Branch("pi0_prob[3]", pi0_prob, "pi0_prob[3]/D");
    tree->Branch("pi0_iteration[3]", pi0_iteration, "pi0_iteration[3]/D");
    tree->Branch("pi0_status[3]", pi0_status, "pi0_status[3]/D");

    tree->Branch("chi2_ref", &chi2_ref);
    tree->Branch("prob_ref", &prob_ref);
    tree->Branch("iteration_ref", &iteration_ref);
    tree->Branch("status_ref", &status_ref);

    tree->Branch("chi2_sig", &chi2_sig);
    tree->Branch("prob_sig", &prob_sig);
    tree->Branch("iteration_sig", &iteration_sig);
    tree->Branch("status_sig", &status_sig);

    tree->Branch("type", &type);
    tree->Branch("truetype", &truetype);

    tree->Branch("decayString",&decayString);

    tree->Branch("kf_beamE",  &kinfittedSig.beamE);
    tree->Branch("kf_gammas", &kinfittedSig.gammas);
    tree->Branch("kf_6g",     &kinfittedSig.etaprimeCand);
    tree->Branch("kf_p",      &kinfittedSig.p);
}

void Etap3pi0::branches::FillKinFit(double beamE, const ant::TParticleList& photons, const ant::TParticlePtr& proton)
{
    assert(photons.size() == 6);

    kinfittedSig.beamE = beamE;
    kinfittedSig.etaprimeCand = TLorentzVector(0,0,0,0);

    for (size_t i = 0 ; i < 6 ; ++i)
    {
        kinfittedSig.gammas.at(i) = *photons.at(i);
        kinfittedSig.etaprimeCand += *photons.at(i);
    }

    kinfittedSig.p = *proton;
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
