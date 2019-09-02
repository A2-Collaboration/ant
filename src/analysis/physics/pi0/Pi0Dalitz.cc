#include "Pi0Dalitz.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "base/WrapTFile.h"

#include <iostream>
#include <memory>

/**
 * Analysis of the pi0->e+e-g channel
 *
 **/

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

Pi0Dalitz::Pi0Dalitz(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get())
{
    CreateHistos();
    AnalysTree.CreateBranches(HistFac.makeTTree("analysis_variables"));
}


void Pi0Dalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    //-- Check the decay string for signal MC pattern
    bool MCpi0eeg  = false;
    string decay;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC))
            decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
    else    decay = "data" + ExpConfig::Setup::Get().GetName();
    if(decay == "(#gamma p) #rightarrow #pi^{0} [ e^{-} #gamma e^{+} ] p ") MCpi0eeg = true;

    //-- MC true stuff
    //--- Particularly for the pi0eeg channel
    TLorentzVector TrueVecpi0;
    TLorentzVector TrueVecep;
    TLorentzVector TrueVecem;
    vector<TLorentzVector> TrueVecGammas;
    if(MCpi0eeg){
        TrueVecpi0 = *utils::ParticleTools::FindParticle(ParticleTypeDatabase::Pi0,event.MCTrue().ParticleTree);
        TrueVecep = *utils::ParticleTools::FindParticle(ParticleTypeDatabase::ePlus,event.MCTrue().ParticleTree);
        TrueVecem = *utils::ParticleTools::FindParticle(ParticleTypeDatabase::eMinus,event.MCTrue().ParticleTree);
        for (const auto& g: utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon,event.MCTrue().ParticleTree)){
            TrueVecGammas.push_back(*g);
        }
    }

    //-- Get full list of candidates, create neutral/charged list
    //   All their info is stored in the tree (add more raw info?)
    const auto& candidates = event.Reconstructed().Candidates;
    TCandidatePtrList recocandidates;
    vector<double> CanCaloE;
    vector<double> CanTheta;
    vector<double> CanPhi;
    vector<double> CanTime;
    vector<double> CanCluSize;
    vector<bool> CanInCB;
    for(const auto& cand : candidates.get_iter()) {
        recocandidates.emplace_back(cand);
        CanCaloE.push_back(cand->CaloEnergy);
        CanTheta.push_back(cand->Theta);
        CanPhi.push_back(cand->Phi);
        CanTime.push_back(cand->Time);
        CanCluSize.push_back(cand->ClusterSize);
        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
            CanInCB.push_back(true);
        else
            CanInCB.push_back(false);
    }

    //-- Loop over the tagger hits and vectors with information related to the tagger
    vector<double> taggweights;
    vector<TLorentzVector> InitialPhotonVec;
    for (const auto &tc : event.Reconstructed().TaggerHits) {
        double cortagtime = triggersimu.GetCorrectedTaggerTime(tc);
        promptrandom.SetTaggerTime(cortagtime);
        double taggweight = promptrandom.FillWeight();
        taggweights.push_back(taggweight);
        InitialPhotonVec.push_back(tc.GetPhotonBeam());
    }

    //-- Fill the tree
    //---- Trigger
    AnalysTree.TBHasTrigg = triggersimu.HasTriggered();
    //---- MC True
    AnalysTree.TBMCpi0eeg = MCpi0eeg;
    AnalysTree.TBTrueVecpi0 = TrueVecpi0;
    AnalysTree.TBTrueVecep = TrueVecep;
    AnalysTree.TBTrueVecem = TrueVecem;
    AnalysTree.TBTrueVecGammas = TrueVecGammas;

    //---- All candidate info
    AnalysTree.TBCanCaloE = CanCaloE;
    AnalysTree.TBCanTheta =  CanTheta;
    AnalysTree.TBCanPhi = CanPhi;
    AnalysTree.TBCanTime = CanTime;
    AnalysTree.TBCanCluSize = CanCluSize;
    AnalysTree.TBCanInCB = CanInCB;

    //--- Tagger info
    AnalysTree.TBPrRndWeight = taggweights;
    AnalysTree.TBInitPhotVec = InitialPhotonVec;

    AnalysTree.Tree->Fill();


}

void Pi0Dalitz::Finish()
{
}

void Pi0Dalitz::ShowResult()
{
/*
    canvas(GetName())
             <<  hTrueGammaE
             << hTrueGammaIM
             << endc; // actually draws the canvas

    ant::canvas(GetName()+": Analysis cuts")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_signal","promptrandom*isSignal")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_background","promptrandom*(!isSignal)")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom","signalcount","0 = bkg, 1 = sig","",BinSettings(2))
            << endc; // actually draws the canvas
*/
}

void Pi0Dalitz::CreateHistos()
{


}
AUTO_REGISTER_PHYSICS(Pi0Dalitz)

