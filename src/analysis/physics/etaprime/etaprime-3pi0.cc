#include "etaprime-3pi0.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "base/std_ext/math.h"
#include "utils/particle_tools.h"

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
    fitter(GetName(), 6)
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
    AddHist1D(cat, "evcount"              , "events after steps", "", "# events", BinSettings(5));

//    fitter.SetupBranches(tree, "EPB");
    fitter.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    signal_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_3Pi0_6g);
    reference_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g);
    bkg_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g);

    vars.SetBranches(tree);

}

TLorentzVector Etap3pi0::MakeLoretzSum(const TParticleList& particles)
{
    TLorentzVector lorenzTemp(0,0,0,0);
    for (const auto& prat: particles)
        lorenzTemp+=*(prat);
    return lorenzTemp;
}








bool Etap3pi0::MakeMCProton(const TEventData& mcdata, TParticlePtr& proton)
{
   const auto& protonlist = mcdata.Particles.Get(ParticleTypeDatabase::Proton);
   if (protonlist.size() != 1)
       return false;
   proton = protonlist.at(0);
   return true;
}

void Etap3pi0::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data   = *event.Reconstructed;
    const auto& mcdata = *event.MCTrue;

/*
    if ( mcdata.ParticleTree )
    {
        if (mcdata.ParticleTree->IsEqual(signal_tree, utils::ParticleTools::MatchByParticleName))
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree).c_str(),1);
        if (mcdata.ParticleTree->IsEqual(reference_tree, utils::ParticleTools::MatchByParticleName))
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree).c_str(),1);
        if (mcdata.ParticleTree->IsEqual(bkg_tree, utils::ParticleTools::MatchByParticleName))
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree).c_str(),1);
    }
*/
    
    const auto& photons            = data.Particles.Get(ParticleTypeDatabase::Photon);
    const auto& protonCandidates   = data.Particles.Get(ParticleTypeDatabase::Proton);
    const auto& mcprotons          = mcdata.Particles.Get(ParticleTypeDatabase::Proton);

    if ( data.Candidates.size() != 7)
        return;

    if (photons.size() != 6)
        return;

    if (protonCandidates.size() != 1)
        return;

    hists.at("steps").at("evcount")->Fill("nParticles",1);

    // geo-cuts

    //proton-cuts -> TAPS

    vars.etaprime = MakeLoretzSum(photons);

    vars.trueProton = *mcprotons.at(0);

    double CBAvgTime = event.Reconstructed->Trigger.CBTiming;
    if(!isfinite(CBAvgTime))
        return;

    for(const TTaggerHit& t : data.TaggerHits )
    {
        promptrandom.SetTaggerHit(t.Time - CBAvgTime);

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        vars.taggWeight    = promptrandom.FillWeight();
        vars.taggE         = t.PhotonEnergy;
        vars.taggCh        = unsigned(t.Channel);
        vars.taggTime      = t.Time;
    }
    tree->Fill();
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
    tree->Branch("proton", &proton);

    tree->Branch("fittedProton",&fittedProton);

    tree->Branch("trueProton", &trueProton);

    tree->Branch("MM", &MM);

    tree->Branch("coplanarity", &coplanarity);

    tree->Branch("taggWeight", &taggWeight);
    tree->Branch("taggE", &taggE);
    tree->Branch("taggCh", &taggCh);
    tree->Branch("taggTime", &taggTime);

    tree->Branch("pi0[3]", &pi0);
    tree->Branch("pi0_chi2[3]", pi0_chi2, "pi0_chi2[3]/D");
    tree->Branch("pi0_prob[3]", pi0_prob, "pi0_prob[3]/D");
    tree->Branch("pi0_iteration[3]", pi0_iteration, "pi0_iteration[3]/D");
    tree->Branch("pi0_status[3]", pi0_status, "pi0_status[3]/D");

    tree->Branch("etaprime", &etaprime);
    tree->Branch("event_chi2", &event_chi2);
    tree->Branch("event_prob", &event_prob);
    tree->Branch("event_iteration", &event_iteration);
    tree->Branch("event_status", &event_status);

    tree->Branch("type", &type);
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
