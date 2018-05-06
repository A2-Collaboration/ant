#include "gp_pi0p.h"

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
 * pi0 photoproduction
 * - decay to 2gamma
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_lheijken_gppi0p::scratch_lheijken_gppi0p(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get())
{
    CreateHistos();
    AnalysTree.CreateBranches(HistFac.makeTTree("analysis_variables"));
}


void scratch_lheijken_gppi0p::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);
    if(!triggersimu.HasTriggered())
        return;

    //-- Check the decay string for signal MC pattern
    bool MCpi0gg  = false;
    string decay;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC))
            decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
    else    decay = "data" + ExpConfig::Setup::Get().GetName();
    if(decay == "(#gamma p) #rightarrow #pi^{0} [ #gamma #gamma ] p ") MCpi0gg = true;

    //-- MC true stuff
    vector<TLorentzVector> GammasTrue;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC)){
        //--- Fetches a list of all gammas in the MCTrue tree
        for (const auto& g: utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon,event.MCTrue().ParticleTree)){
            GammasTrue.push_back(*g);
        }
    }
    //--- Particularly for the pi0gg channel
    TLorentzVector TrueVecpi0;
    if(MCpi0gg){
        TrueVecpi0 = *utils::ParticleTools::FindParticle(ParticleTypeDatabase::Pi0,event.MCTrue().ParticleTree);
    }

    //-- Get full list of particles, create neutral/charged list
    //   All their info is stored in the tree (add more raw info?)
    const auto& candidates = event.Reconstructed().Candidates;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    vector<double> NeuCanCaloE;
    vector<double> NeuCanTheta;
    vector<double> NeuCanPhi;
    vector<double> NeuCanTime;
    vector<double> NeuCanCluSize;
    vector<bool> NeuCanInCB;
    vector<double> ChaCanCaloE;
    vector<double> ChaCanVetoE;
    vector<double> ChaCanTheta;
    vector<double> ChaCanPhi;
    vector<double> ChaCanTime;
    vector<double> ChaCanCluSize;
    vector<bool> ChaCanInCB;
    vector<bool> ChaCanInPID;
    for(const auto& cand : candidates.get_iter()) {
        if(cand->VetoEnergy == 0.0) {
            neutral.emplace_back(cand);
            NeuCanCaloE.push_back(cand->CaloEnergy);
            NeuCanTheta.push_back(cand->Theta);
            NeuCanPhi.push_back(cand->Phi);
            NeuCanTime.push_back(cand->Time);
            NeuCanCluSize.push_back(cand->ClusterSize);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                NeuCanInCB.push_back(true);
            else
                NeuCanInCB.push_back(false);
        }
        else{
            charged.emplace_back(cand);
            neutral.emplace_back(cand);
            ChaCanCaloE.push_back(cand->CaloEnergy);
            ChaCanVetoE.push_back(cand->VetoEnergy);
            ChaCanTheta.push_back(cand->Theta);
            ChaCanPhi.push_back(cand->Phi);
            ChaCanTime.push_back(cand->Time);
            ChaCanCluSize.push_back(cand->ClusterSize);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                ChaCanInCB.push_back(true);
            else
                ChaCanInCB.push_back(false);
            if(cand->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID)
                ChaCanInPID.push_back(true);
            else
                ChaCanInPID.push_back(false);
        }
    }

    //-- Focus on pi->gg events
    //--- the variables to be saved to the tree
    vector<double> pi0gg_times;
    vector<bool> pi0ggOnlyCB;
    vector<TLorentzVector> pi0ggVecRec;
    vector<int> pi0ggNeuCombInd1;
    vector<int> pi0ggNeuCombInd2;
    //--- the variables used to fill the histograms
    double pi0gg_IM=-100;
    vector<double> pi0ggEg;
    vector<double> pi0ggTheg;
    vector<double> pi0ggPhig;
    double pi0ggEpi0 = -100;
    double pi0ggThepi0 = -100;
    double pi0ggPhipi0 = -100;
    //--- Only 2 neutrals
    if(NeuCanCaloE.size() == 2){
        bool ispi0ggOnlyCB = true;
        TParticle pi0ggOnly2g;
        double pi0ggOnly2g_time = 0;
        for (const auto& photon : neutral) {
            pi0ggOnly2g+=TParticle(ParticleTypeDatabase::Photon, photon);
            pi0ggOnly2g_time+=photon->Time;
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                ispi0ggOnlyCB = false;
            pi0ggEg.push_back(photon->CaloEnergy);
            pi0ggTheg.push_back(photon->Theta);
            pi0ggPhig.push_back(photon->Phi);
        }
        pi0ggOnly2g_time = pi0ggOnly2g_time/2.;
        pi0gg_times.push_back(pi0ggOnly2g_time);
        pi0ggOnlyCB.push_back(ispi0ggOnlyCB);
        pi0ggNeuCombInd1.push_back(0);
        pi0ggNeuCombInd2.push_back(1);
        pi0ggVecRec.push_back(pi0ggOnly2g);
        pi0gg_IM = pi0ggOnly2g.M();
        pi0ggEpi0 = pi0ggOnly2g.E;
        pi0ggThepi0 = pi0ggOnly2g.Theta();
        pi0ggPhipi0 = pi0ggOnly2g.Phi();
    }

    //-- Focus on pi->e+e-g events
    //--- the variables to be saved to the tree
    vector<double> pi0DD_times;
    vector<bool> pi0DDOnlyCB;
    vector<TLorentzVector> pi0DDVecRec;
    //--- the variables used to fill the histograms
    double pi0DD_IM=-100;
    double pi0DDEpi0 = -100;
    double pi0DDThepi0 = -100;
    double pi0DDPhipi0 = -100;
    //--- exactly 1 neutral and 2 charged
    if(NeuCanCaloE.size() == 1 && ChaCanCaloE.size() == 2){
        bool ispi0DDOnlyCB = true;
        TParticle pi0DD1Ne2Ch;
        double pi0DD1Ne2Ch_time = 0;
        for (const auto& photon : neutral) { // there should be only one though, couldn't think of a better way to extract the first (and only) neutral
            pi0DD1Ne2Ch+=TParticle(ParticleTypeDatabase::Photon, photon);
            pi0DD1Ne2Ch_time+=photon->Time;
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                ispi0DDOnlyCB = false;
        }
        for (const auto& lepton : charged) {
            pi0DD1Ne2Ch+=TParticle(ParticleTypeDatabase::eMinus, lepton);
            pi0DD1Ne2Ch_time+=lepton->Time;
            if(lepton->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                ispi0DDOnlyCB = false;
        }
        pi0DD1Ne2Ch_time = pi0DD1Ne2Ch_time/3.;
        pi0DD_times.push_back(pi0DD1Ne2Ch_time);
        pi0DDOnlyCB.push_back(ispi0DDOnlyCB);
        pi0DDVecRec.push_back(pi0DD1Ne2Ch);
        pi0DD_IM = pi0DD1Ne2Ch.M();
        pi0DDEpi0 = pi0DD1Ne2Ch.E;
        pi0DDThepi0 = pi0DD1Ne2Ch.Theta();
        pi0DDPhipi0 = pi0DD1Ne2Ch.Phi();
    }


    //-- Loop over the tagger hits and fill the histograms and the tree
    for (const auto &tc : event.Reconstructed().TaggerHits) {
        double cortagtime = triggersimu.GetCorrectedTaggerTime(tc);
        promptrandom.SetTaggerTime(cortagtime);
        double taggweight = promptrandom.FillWeight();
        TLorentzVector InitialPhotonVec = tc.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());

        //--- Histograms
        //---- Checks
        hPromRandWei->Fill(taggweight);
        hTaggTimeChannel->Fill(tc.Time,tc.Channel);
        hTaggCorTimeChannel->Fill(cortagtime,tc.Channel);
        //---- MC True
        for(const auto& tg : GammasTrue){
            hTrueGammaE->Fill(tg.E(),taggweight);
            hTrueThevsEg->Fill(tg.Theta()*radtodeg,tg.E(),taggweight);
            hTrueThevsPhig->Fill(tg.Theta()*radtodeg,tg.Phi()*radtodeg,taggweight);
        }
        if(MCpi0gg){
            hTrueIMgg->Fill((GammasTrue[0]+GammasTrue[1]).M(),taggweight);
            hTrueMMgg->Fill((InitialPhotonVec + InitialProtonVec - GammasTrue[0] - GammasTrue[1]).M(), taggweight);
            hTrueThevsEpi0->Fill(TrueVecpi0.Theta()*radtodeg,TrueVecpi0.E(),taggweight);
            hTrueThevsPhipi0->Fill(TrueVecpi0.Theta()*radtodeg,TrueVecpi0.Phi()*radtodeg,taggweight);
        }
        //---- All candidate info
        for(unsigned int nci=0; nci<NeuCanTheta.size(); nci++){
            hNeuCanThevsCaloE->Fill(NeuCanTheta[nci]*radtodeg,NeuCanCaloE[nci],taggweight);
            hNeuCanThevsPhi->Fill(NeuCanTheta[nci]*radtodeg,NeuCanPhi[nci]*radtodeg,taggweight);
            if(NeuCanInCB[nci]){
                hNeuCanCBTime->Fill(NeuCanTime[nci],taggweight);
                hNeuCanCBCluSize->Fill(NeuCanCluSize[nci],taggweight);
            }
            else{
                hNeuCanTAPSTime->Fill(NeuCanTime[nci],taggweight);
                hNeuCanTAPSCluSize->Fill(NeuCanCluSize[nci],taggweight);
            }
        }
        for(unsigned int cci=0; cci<ChaCanTheta.size(); cci++){
            hChaCanThevsCaloE->Fill(ChaCanTheta[cci]*radtodeg,ChaCanCaloE[cci],taggweight);
            hChaCanThevsVetoE->Fill(ChaCanTheta[cci]*radtodeg,ChaCanVetoE[cci],taggweight);
            hChaCanThevsPhi->Fill(ChaCanTheta[cci]*radtodeg,ChaCanPhi[cci]*radtodeg,taggweight);
            if(ChaCanInCB[cci]){
                hChaCanCBTime->Fill(ChaCanTime[cci],taggweight);
                hChaCanCBCluSize->Fill(ChaCanCluSize[cci],taggweight);
            }
            else{
                hChaCanTAPSTime->Fill(ChaCanTime[cci],taggweight);
                hChaCanTAPSCluSize->Fill(ChaCanCluSize[cci],taggweight);
            }
        }
        //---- pi0gg stuff
        if(NeuCanCaloE.size() == 2){
            if(pi0ggOnlyCB[0]){
                hpi0ggO2gOCB_time->Fill(pi0gg_times[0],taggweight);
                hpi0ggO2gOCB_IM->Fill(pi0gg_IM,taggweight);
                hpi0ggO2gOCB_IMggVsTagCh->Fill(pi0gg_IM,tc.Channel,taggweight);
                hpi0ggO2gOCB_TagTimeVsChan->Fill(tc.Time,tc.Channel);
                hpi0ggO2gOCB_TagCorTimeVsChan->Fill(cortagtime,tc.Channel);
            }
            else {
                hpi0ggO2gCBTA_time->Fill(pi0gg_times[0],taggweight);
                hpi0ggO2gCBTA_IM->Fill(pi0gg_IM,taggweight);
                hpi0ggO2gCBTA_IMggVsTagCh->Fill(pi0gg_IM,tc.Channel,taggweight);
                hpi0ggO2gCBTA_TagTimeVsChan->Fill(tc.Time,tc.Channel);
                hpi0ggO2gCBTA_TagCorTimeVsChan->Fill(cortagtime,tc.Channel);
            }
            for(int it=0; it<2; it++){
                hpi0ggO2g_gThevsE->Fill(pi0ggTheg[it]*radtodeg,pi0ggEg[it],taggweight);
                hpi0ggO2g_gThevsPhi->Fill(pi0ggTheg[it]*radtodeg,pi0ggPhig[it]*radtodeg,taggweight);
            }
            hpi0ggO2g_pi0ThevsE->Fill(pi0ggThepi0*radtodeg,pi0ggEpi0,taggweight);
            hpi0ggO2g_pi0ThevsPhi->Fill(pi0ggThepi0*radtodeg,pi0ggPhipi0*radtodeg,taggweight);
        }
        //---- pi0DD stuff
        if(NeuCanCaloE.size() == 1 && ChaCanCaloE.size() == 2){
            if(pi0DDOnlyCB[0]){
                hpi0DD1N2COCB_time->Fill(pi0DD_times[0],taggweight);
                hpi0DD1N2COCB_IM->Fill(pi0DD_IM,taggweight);
            }
            else {
                hpi0DD1N2CCBTA_time->Fill(pi0DD_times[0],taggweight);
                hpi0DD1N2CCBTA_IM->Fill(pi0DD_IM,taggweight);
            }
            hpi0DD1N2C_pi0ThevsE->Fill(pi0DDThepi0*radtodeg,pi0DDEpi0,taggweight);
            hpi0DD1N2C_pi0ThevsPhi->Fill(pi0DDThepi0*radtodeg,pi0DDPhipi0*radtodeg,taggweight);
        }



        //--- Tree
        //---- Checks
        AnalysTree.TBPrRndWeight = promptrandom.FillWeight();
        //---- MC True
        AnalysTree.TBMCpi0gg = MCpi0gg;
        AnalysTree.TBTrueGammas = GammasTrue;
        AnalysTree.TBTrueVecpi0 = TrueVecpi0;
        //---- All candidate info
        AnalysTree.TBNeuCanCaloE = NeuCanCaloE;
        AnalysTree.TBNeuCanTheta =  NeuCanTheta;
        AnalysTree.TBNeuCanPhi = NeuCanPhi;
        AnalysTree.TBNeuCanTime = NeuCanTime;
        AnalysTree.TBNeuCanCluSize = NeuCanCluSize;
        AnalysTree.TBNeuCanInCB = NeuCanInCB;
        AnalysTree.TBChaCanCaloE = ChaCanCaloE;
        AnalysTree.TBChaCanVetoE = ChaCanVetoE;
        AnalysTree.TBChaCanTheta = ChaCanTheta;
        AnalysTree.TBChaCanPhi = ChaCanPhi;
        AnalysTree.TBChaCanTime = ChaCanTime;
        AnalysTree.TBChaCanCluSize = ChaCanCluSize;
        AnalysTree.TBChaCanInCB = ChaCanInCB;
        AnalysTree.TBChaCanInPID = ChaCanInPID;
        //---- pi0gg stuff
        AnalysTree.TBpi0ggTime = pi0gg_times;
        AnalysTree.TBpi0ggVecRec = pi0ggVecRec;
        AnalysTree.TBpi0ggOnlyCB = pi0ggOnlyCB;
        AnalysTree.TBpi0ggNeuCombInd1 = pi0ggNeuCombInd1;
        AnalysTree.TBpi0ggNeuCombInd2 = pi0ggNeuCombInd2;
        //---- pi0DD stuff
        AnalysTree.TBpi0DDTime = pi0DD_times;
        AnalysTree.TBpi0DDVecRec = pi0DDVecRec;
        AnalysTree.TBpi0DDOnlyCB = pi0DDOnlyCB;

        AnalysTree.Tree->Fill();

    }

}

void scratch_lheijken_gppi0p::Finish()
{
}

void scratch_lheijken_gppi0p::ShowResult()
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

void scratch_lheijken_gppi0p::CreateHistos()
{
    auto hfTrueMC = new HistogramFactory("TrueMC", HistFac, "");
    auto hfChecksH = new HistogramFactory("ChecksHistos", HistFac, "");
    auto hfAllCandH = new HistogramFactory("AllCandHistos", HistFac, "");
    auto hfPi0ggH = new HistogramFactory("Pi0ggHistos", HistFac, "");
    auto hfPi0DDH = new HistogramFactory("Pi0DDHistos", HistFac, "");


    const BinSettings ETGBins(100,0.,500.);
    const BinSettings PhiBins(500,-180.,180.);
    const BinSettings TheBins(500,0.,180.);
    const BinSettings TimeBins(1000,-199.5,200.5);
    const BinSettings TimeBins1ns(400,-200,200);
    const BinSettings EnergyVetoBins(100,-0.5,49.5);
    const BinSettings EnergySBins(500,-0.5,499.5);
    const BinSettings EnergyMBins(1000,-0.5,999.5);
    const BinSettings EnergyLBins(4000,-0.5,3999.5);
    const BinSettings IMggBins(100,0.,500.);
    const BinSettings MMpi0Bins(1000,-0.5,1999.5);

    //--- MC True
    hTrueGammaE = hfTrueMC->makeTH1D("energies of true gammas","E_{#gamma}","",ETGBins,"TrueGammaE", true);
    hTrueIMgg = hfTrueMC->makeTH1D("True IMgg","IM(#gamma#gamma)","",IMggBins,"TrueIMgg",true);
    hTrueMMgg = hfTrueMC->makeTH1D("True MMpp","MM(#gamma#gamma)","",MMpi0Bins,"TrueMMgg",true);
    hTrueThevsEg = hfTrueMC->makeTH2D("True Theta vs E for photons","#theta_{#gamma}","E_{#gamma}",TheBins, EnergyMBins,"TrueThevsEg",true);
    hTrueThevsPhig = hfTrueMC->makeTH2D("True Theta vs Phi for photons","#theta_{#gamma}","#phi_{#gamma}",TheBins, PhiBins,"TrueThevsPhig",true);
    hTrueThevsEpi0 = hfTrueMC->makeTH2D("True Theta vs E for pi0","#theta_{#pi^{0}}","E_{#pi^{0}}",TheBins, EnergyMBins,"TrueThevsEpi0",true);
    hTrueThevsPhipi0 = hfTrueMC->makeTH2D("True Theta vs Phi for pi0","#theta_{#pi^{0}}","#phi_{#pi^{0}}",TheBins, PhiBins,"TrueThevsPhipi0",true);

    //--- Checks
    hPromRandWei = hfChecksH->makeTH1D("Promt Rnd weights","w","",{401,-2,2},"PromRandWei",true);
    hTaggTimeChannel = hfChecksH->makeTH2D("Tagger time vs channel","Tagger time","Tagger channel",TimeBins1ns,BinSettings(352),"hTaggTimeChannel",true);
    hTaggCorTimeChannel = hfChecksH->makeTH2D("Tagger corr. time vs channel","Tagger time - CBsum time","Tagger channel",TimeBins1ns,BinSettings(352),"hTaggCorTimeChannel",true);

    //--- All candidate info
    hNeuCanThevsCaloE = hfAllCandH->makeTH2D("NeuCand Theta vs CaloE","#theta","E",TheBins, EnergySBins,"NeuCanThevsCaloE", true);
    hNeuCanThevsPhi = hfAllCandH->makeTH2D("NeuCand Theta vs Phi","#theta","#phi",TheBins, PhiBins,"NeuCanThevsPhi", true);
    hNeuCanCBTime = hfAllCandH->makeTH1D("NeuCand CB Time","Time","",TimeBins,"NeuCanCBTime", true);
    hNeuCanTAPSTime = hfAllCandH->makeTH1D("NeuCand TAPS Time","Time","",TimeBins,"NeuCanTAPSTime", true);
    hNeuCanCBCluSize = hfAllCandH->makeTH1D("NeuCand CB CluSize","No. crystals","",{30,-0.5,29.5},"NeuCanCBCluSize",true);
    hNeuCanTAPSCluSize = hfAllCandH->makeTH1D("NeuCand TAPS CluSize","No. crystals","",{30,-0.5,29.5},"NeuCanTAPSCluSize",true);
    hChaCanThevsCaloE = hfAllCandH->makeTH2D("ChaCand Theta vs CaloE","#theta","E",TheBins, EnergySBins,"ChaCanThevsCaloE", true);
    hChaCanThevsVetoE = hfAllCandH->makeTH2D("ChaCand Theta vs VetoE","#theta","E",TheBins, EnergyVetoBins,"ChaCanThevsVetoE", true);
    hChaCanThevsPhi = hfAllCandH->makeTH2D("ChaCand Theta vs Phi","#theta","#phi",TheBins, PhiBins,"ChaCanThevsPhi", true);
    hChaCanCBTime = hfAllCandH->makeTH1D("ChaCand CB Time","Time","",TimeBins,"ChaCanCBTime", true);
    hChaCanTAPSTime = hfAllCandH->makeTH1D("ChaCand TAPS Time","Time","",TimeBins,"ChaCanTAPSTime", true);
    hChaCanCBCluSize = hfAllCandH->makeTH1D("ChaCand CB CluSize","No. crystals","",{30,-0.5,29.5},"ChaCanCBCluSize",true);
    hChaCanTAPSCluSize = hfAllCandH->makeTH1D("ChaCand TAPS CluSize","No. crystals","",{30,-0.5,29.5},"ChaCanTAPSCluSize",true);

    //--- pi0gg stuff
    hpi0ggO2gOCB_time = hfPi0ggH->makeTH1D("pi0gg av.Time, only 2g, only CB","T","",TimeBins,"pi0ggO2gOCB_time", true);
    hpi0ggO2gCBTA_time = hfPi0ggH->makeTH1D("pi0gg av.Time, only 2g, CB and TAPS","T","",TimeBins,"pi0ggO2gCBTA_time", true);
    hpi0ggO2gOCB_IM = hfPi0ggH->makeTH1D("pi0gg IMgg, only 2, only CB","IM(#gamma#gamma)","",IMggBins,"pi0ggO2gOCB_IM",true);
    hpi0ggO2gCBTA_IM = hfPi0ggH->makeTH1D("pi0gg IMgg, only 2, CB and TAPS","IM(#gamma#gamma)","",IMggBins,"pi0ggO2gCBTA_IM",true);
    hpi0ggO2g_gThevsE = hfPi0ggH->makeTH2D("pi0gg Theta vs E for photons, only 2g","#theta_{#gamma}","E_{#gamma}",TheBins,EnergySBins,"pi0ggO2g_gThevsE",true);
    hpi0ggO2g_gThevsPhi = hfPi0ggH->makeTH2D("pi0gg Theta vs Phi for photons, only 2g","#theta_{#gamma}","#phi_{#gamma}",TheBins,PhiBins,"pi0ggO2g_gThevsPhi",true);
    hpi0ggO2g_pi0ThevsE = hfPi0ggH->makeTH2D("pi0gg Theta vs E for pi0, only 2g","#theta_{#pi^{0}}","E_{#pi^{0}}",TheBins,EnergyMBins,"pi0ggO2g_pi0ThevsE",true);
    hpi0ggO2g_pi0ThevsPhi = hfPi0ggH->makeTH2D("pi0gg Theta vs Phi for pi0, only 2g","#theta_{#pi^{0}}","#phi_{#pi^{0}}",TheBins,PhiBins,"pi0ggO2g_pi0ThevsPhi",true);
    hpi0ggO2gOCB_IMggVsTagCh = hfPi0ggH->makeTH2D("pi0gg IMgg vs Tagger channel, only 2g, only CB","IM(#gamma#gamma)","Tagger channel",IMggBins,BinSettings(352),"pi0ggO2gOCB_IMggVsTagCh",true);
    hpi0ggO2gCBTA_IMggVsTagCh = hfPi0ggH->makeTH2D("pi0gg IMgg vs Tagger channel, only 2g, CB and TAPS","IM(#gamma#gamma)","Tagger channel",IMggBins,BinSettings(352),"pi0ggO2gCBTA_IMggVsTagCh",true);
    hpi0ggO2gOCB_TagTimeVsChan = hfPi0ggH->makeTH2D("pi0gg Tagger time vs channel, only 2g, only CB","Tagger time","Tagger channel",TimeBins1ns,BinSettings(352),"hpi0ggO2gOCB_TagTimeVsChan",true);
    hpi0ggO2gCBTA_TagTimeVsChan = hfPi0ggH->makeTH2D("pi0gg Tagger time vs channel, only 2g, CB and TAPS","Tagger time","Tagger channel",TimeBins1ns,BinSettings(352),"hpi0ggO2gCBTA_TagTimeVsChan",true);
    hpi0ggO2gOCB_TagCorTimeVsChan = hfPi0ggH->makeTH2D("pi0gg Tagger corr. time vs channel, only 2g, only CB","Tagger time - CBsum time","Tagger channel",TimeBins1ns,BinSettings(352),"hpi0ggO2gOCB_TagCorTimeVsChan",true);
    hpi0ggO2gCBTA_TagCorTimeVsChan = hfPi0ggH->makeTH2D("pi0gg Tagger corr. time vs channel, only 2g, CB and TAPS","Tagger time - CBsum time","Tagger channel",TimeBins1ns,BinSettings(352),"hpi0ggO2gCBTA_TagCorTimeVsChan",true);


    //---pi0DD stuff
    hpi0DD1N2COCB_time = hfPi0DDH->makeTH1D("pi0e+e-g av.Time, only 1Ne2Ch, only CB","T","",TimeBins,"pi0DD1N2COCB_time", true);
    hpi0DD1N2CCBTA_time = hfPi0DDH->makeTH1D("pi0e+e-g av.Time, only 1Ne2Ch, CB and TAPS","T","",TimeBins,"pi0DD1N2CCBTA_time", true);
    hpi0DD1N2COCB_IM = hfPi0DDH->makeTH1D("pi0e+e-g IMeeg, only 1Ne2Ch, only CB","IM(e^{+}e^{-}#gamma)","",IMggBins,"pi0DD1N2COCB_IM",true);
    hpi0DD1N2CCBTA_IM = hfPi0DDH->makeTH1D("pi0e+e-g IMeeg, only 1Ne2Ch, CB and TAPS","IM(e^{+}e^{-}#gamma)","",IMggBins,"pi0DD1N2CCBTA_IM",true);
    hpi0DD1N2C_pi0ThevsE = hfPi0DDH->makeTH2D("pi0e+e-g Theta vs E for pi0, only 1Ne2Ch","#theta_{#pi^{0}}","E_{#pi^{0}}",TheBins,EnergyMBins,"pi0DD1N2C_pi0ThevsE",true);
    hpi0DD1N2C_pi0ThevsPhi = hfPi0DDH->makeTH2D("pi0e+e-g Theta vs Phi for pi0, only 1Ne2Ch","#theta_{#pi^{0}}","#phi_{#pi^{0}}",TheBins,PhiBins,"pi0DD1N2C_pi0ThevsPhi",true);
}
AUTO_REGISTER_PHYSICS(scratch_lheijken_gppi0p)

