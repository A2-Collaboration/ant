#include "gp_wp_pi0g.h"

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
#include "base/WrapTFile.h"

#include <iostream>
#include <memory>

/**
 * A gp->pw->pi0g class (only a first, very simplistic version)
**/

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_lheijken_gpwppi0g::scratch_lheijken_gpwppi0g(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get())
{
    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    nTagger = tagger_detector->GetNChannels();

    CreateHistos();
}


void scratch_lheijken_gpwppi0g::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);
    if(!triggersimu.HasTriggered())
        return;

    //-- Check the decay string for signal MC pattern
    bool MCwpi0g  = false;
    string decay;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC)){
        const auto& particletree = event.MCTrue().ParticleTree;
        if (!particletree) {
            LOG(ERROR) << "No particle tree found, only Geant or Pluto file provided, not both";
            return;
        }
        // check if the current event is the signal
        MCwpi0g = particletree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g),
                                                    utils::ParticleTools::MatchByParticleName);
        decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
    }
    else
        decay = "data" + ExpConfig::Setup::Get().GetName();

    //-- MC true stuff
    vector<TLorentzVector> GammasTrue;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC)){
        //--- Fetches a list of all gammas in the MCTrue tree
        for (const auto& g: utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon,event.MCTrue().ParticleTree)){
            GammasTrue.push_back(*g);
        }
    }
    //--- Particularly for the wpi0g channel
    TLorentzVector TrueVecw;
    if(MCwpi0g){
        TrueVecw = *utils::ParticleTools::FindParticle(ParticleTypeDatabase::Omega,event.MCTrue().ParticleTree);
    }


    //-- Get full list of particles, create neutral/charged list
    str_recocand recocand;
    const auto& candidates = event.Reconstructed().Candidates;
    for(const auto& cand : candidates.get_iter()) {
        if(cand->VetoEnergy == 0.0) {
            recocand.neutral.emplace_back(cand);
            recocand.NeuCanCaloE.push_back(cand->CaloEnergy);
            recocand.NeuCanTheta.push_back(cand->Theta);
            recocand.NeuCanPhi.push_back(cand->Phi);
            recocand.NeuCanTime.push_back(cand->Time);
            recocand.NeuCanCluSize.push_back(cand->ClusterSize);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                recocand.NeuCanInCB.push_back(true);
            else
                recocand.NeuCanInCB.push_back(false);
        }
        else{
            recocand.charged.emplace_back(cand);
            recocand.ChaCanCaloE.push_back(cand->CaloEnergy);
            recocand.ChaCanVetoE.push_back(cand->VetoEnergy);
            recocand.ChaCanTheta.push_back(cand->Theta);
            recocand.ChaCanPhi.push_back(cand->Phi);
            recocand.ChaCanTime.push_back(cand->Time);
            recocand.ChaCanCluSize.push_back(cand->ClusterSize);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                recocand.ChaCanInCB.push_back(true);
            else
                recocand.ChaCanInCB.push_back(false);
            if(cand->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID)
                recocand.ChaCanInPID.push_back(true);
            else
                recocand.ChaCanInPID.push_back(false);
        }
    }

    //-- Focus on w->pi0g events ---------------------
    //--- the variables used to fill the histograms
    str_wpi0gO3g wpi0gO3gInfo;
    double wpi0gOnly3g_time = -1000;
    bool iswpi0gOnlyCB = true;
    TParticle wpi0gOnly3g;
    //--- Only 3 neutrals
    if(recocand.NeuCanCaloE.size() == 3){
        wpi0gOnly3g_time=0;
        for (const auto& photon : recocand.neutral) {
            wpi0gOnly3g+=TParticle(ParticleTypeDatabase::Photon, photon);
            wpi0gOnly3g_time+=photon->Time;
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                iswpi0gOnlyCB = false;
            TLorentzVector TVg = TParticle(ParticleTypeDatabase::Photon, photon);
            wpi0gO3gInfo.LV3g.push_back(TVg);
        }
        wpi0gO3gInfo.LVw = wpi0gOnly3g;
        wpi0gO3gInfo.avgT3g = wpi0gOnly3g_time/3.;
        wpi0gO3gInfo.onlyinCB = iswpi0gOnlyCB;
        //--- And only 1 charged
        if(recocand.ChaCanCaloE.size() == 1){
            wpi0gO3gInfo.LVp = TParticle(ParticleTypeDatabase::Proton, recocand.charged.at(0));
        }
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
        if(MCwpi0g){
            hTrueIMggg->Fill((GammasTrue[0]+GammasTrue[1]+GammasTrue[2]).M(),taggweight);
            hTrueMMp->Fill((InitialPhotonVec + InitialProtonVec - GammasTrue[0] - GammasTrue[1] - GammasTrue[2]).M(), taggweight);
            hTrueThevsEw->Fill(TrueVecw.Theta()*radtodeg,TrueVecw.E(),taggweight);
            hTrueThevsPhiw->Fill(TrueVecw.Theta()*radtodeg,TrueVecw.Phi()*radtodeg,taggweight);
        }

        //---- All candidate info
        FillRecoCandHistos(recocand,taggweight);

        //---- wpi0g stuff
        if(recocand.NeuCanCaloE.size() == 3){
            FillO3gHistos(wpi0gO3gInfo, InitialPhotonVec, tc, cortagtime, taggweight);

            //----- If there's only one charged candidate (undevelopped part!)
            if(recocand.ChaCanCaloE.size() == 1){
                hO3gO1p_IMggg->Fill((wpi0gO3gInfo.LVw).M(),taggweight);
                hO3gO1p_IMggg_Thggg->Fill((wpi0gO3gInfo.LVw).M(),(wpi0gO3gInfo.LVw).Theta()*radtodeg,taggweight);
                hO3gO1p_MM->Fill((InitialPhotonVec+InitialProtonVec-wpi0gO3gInfo.LVw).M(),taggweight);
            }
        }
    }
}

void scratch_lheijken_gpwppi0g::Finish()
{
}

void scratch_lheijken_gpwppi0g::ShowResult()
{
}

void scratch_lheijken_gpwppi0g::CreateHistos()
{
    auto hfTrueMC = new HistogramFactory("TrueMC", HistFac, "");
    auto hfChecksH = new HistogramFactory("ChecksHistos", HistFac, "");
    auto hfAllCandH = new HistogramFactory("AllCandHistos", HistFac, "");
    auto hfPi0gH = new HistogramFactory("Pi0gHistos", HistFac, "");

    const BinSettings ETGBins(400,0.,1200.);
    const BinSettings PhiBins(500,-190.,190.);
    const BinSettings TheBins(500,0.,190.);
    const BinSettings cosThBins(101,-1.,1.);
    const BinSettings TimeBins(1000,-199.5,200.5);
    const BinSettings TimeBins1ns(400,-200,200);
    const BinSettings EnergyVetoBins(100,-0.5,49.5);
    const BinSettings EnergySBins(800,-0.5,799.5);
    const BinSettings EnergyMBins(1500,-0.5,1499.5);
    const BinSettings EnergyLBins(4000,-0.5,3999.5);
    const BinSettings IMgggBins(200,0.,1000.);
    const BinSettings MMpBins(500,-0.5,1999.5);

    //--- MC True
    hTrueGammaE = hfTrueMC->makeTH1D("energies of true gammas","E_{#gamma}","",ETGBins,"TrueGammaE", true);
    hTrueIMggg = hfTrueMC->makeTH1D("True IMggg","IM(#gamma#gamma#gamma)","",IMgggBins,"TrueIMggg",true);
    hTrueMMp = hfTrueMC->makeTH1D("True MMp","MM(p)","",MMpBins,"TrueMMp",true);
    hTrueThevsEg = hfTrueMC->makeTH2D("True Theta vs E for photons","#theta_{#gamma}","E_{#gamma}",TheBins, EnergyMBins,"TrueThevsEg",true);
    hTrueThevsPhig = hfTrueMC->makeTH2D("True Theta vs Phi for photons","#theta_{#gamma}","#phi_{#gamma}",TheBins, PhiBins,"TrueThevsPhig",true);
    hTrueThevsEw = hfTrueMC->makeTH2D("True Theta vs E for w","#theta_{#omega}}","E_{#omega}",TheBins, EnergyMBins,"TrueThevsEw",true);
    hTrueThevsPhiw = hfTrueMC->makeTH2D("True Theta vs Phi for w","#theta_{#omega}","#phi_{#omega}",TheBins, PhiBins,"TrueThevsPhiw",true);

    //--- Checks
    hPromRandWei = hfChecksH->makeTH1D("Promt Rnd weights","w","",{401,-2,2},"PromRandWei",true);
    hTaggTimeChannel = hfChecksH->makeTH2D("Tagger time vs channel","Tagger time","Tagger channel",TimeBins1ns,BinSettings(nTagger+20),"hTaggTimeChannel",true);
    hTaggCorTimeChannel = hfChecksH->makeTH2D("Tagger corr. time vs channel","Tagger time - CBsum time","Tagger channel",TimeBins1ns,BinSettings(nTagger+20),"hTaggCorTimeChannel",true);

    //--- All candidate info
    hNeuCanThevsCaloE = hfAllCandH->makeTH2D("NeuCand Theta vs CaloE","#theta","E",TheBins, EnergyMBins,"NeuCanThevsCaloE", true);
    hNeuCanThevsPhi = hfAllCandH->makeTH2D("NeuCand Theta vs Phi","#theta","#phi",TheBins, PhiBins,"NeuCanThevsPhi", true);
    hNeuCanCBTime = hfAllCandH->makeTH1D("NeuCand CB Time","Time","",TimeBins,"NeuCanCBTime", true);
    hNeuCanTAPSTime = hfAllCandH->makeTH1D("NeuCand TAPS Time","Time","",TimeBins,"NeuCanTAPSTime", true);
    hNeuCanCBCluSize = hfAllCandH->makeTH1D("NeuCand CB CluSize","No. crystals","",{30,-0.5,29.5},"NeuCanCBCluSize",true);
    hNeuCanTAPSCluSize = hfAllCandH->makeTH1D("NeuCand TAPS CluSize","No. crystals","",{30,-0.5,29.5},"NeuCanTAPSCluSize",true);
    hChaCanThevsCaloE = hfAllCandH->makeTH2D("ChaCand Theta vs CaloE","#theta","E",TheBins, EnergyMBins,"ChaCanThevsCaloE", true);
    hChaCanThevsVetoE = hfAllCandH->makeTH2D("ChaCand Theta vs VetoE","#theta","E",TheBins, EnergyVetoBins,"ChaCanThevsVetoE", true);
    hChaCanThevsPhi = hfAllCandH->makeTH2D("ChaCand Theta vs Phi","#theta","#phi",TheBins, PhiBins,"ChaCanThevsPhi", true);
    hChaCanCBTime = hfAllCandH->makeTH1D("ChaCand CB Time","Time","",TimeBins,"ChaCanCBTime", true);
    hChaCanTAPSTime = hfAllCandH->makeTH1D("ChaCand TAPS Time","Time","",TimeBins,"ChaCanTAPSTime", true);
    hChaCanCBCluSize = hfAllCandH->makeTH1D("ChaCand CB CluSize","No. crystals","",{30,-0.5,29.5},"ChaCanCBCluSize",true);
    hChaCanTAPSCluSize = hfAllCandH->makeTH1D("ChaCand TAPS CluSize","No. crystals","",{30,-0.5,29.5},"ChaCanTAPSCluSize",true);

    //--- w->pi0g stuff
    hO3gOCB_time = hfPi0gH->makeTH1D("wpi0g av.Time, only 3g, only CB","T","",TimeBins,"O3gOCB_time", true);
    hO3gCBTA_time = hfPi0gH->makeTH1D("wpi0g av.Time, only 3g, CB and TAPS","T","",TimeBins,"O3gCBTA_time", true);
    hO3gOCB_IM = hfPi0gH->makeTH1D("wpi0g IMggg, only 3g, only CB","IM(#gamma#gamma#gamma)","",IMgggBins,"O3gOCB_IM",true);
    hO3gCBTA_IM = hfPi0gH->makeTH1D("wpi0g IMggg, only 3g, CB and TAPS","IM(#gamma#gamma#gamma)","",IMgggBins,"O3gCBTA_IM",true);
    hO3gOCB_MM = hfPi0gH->makeTH1D("wpi0g MMggg, only 3g, only CB","MM(#gamma#gamma#gamma)","",MMpBins,"O3gOCB_MM",true);
    hO3gCBTA_MM = hfPi0gH->makeTH1D("wpi0g MMggg, only 3g, CB and TAPS","MM(#gamma#gamma#gamma)","",MMpBins,"O3gCBTA_MM",true);
    hO3g_gThevsE = hfPi0gH->makeTH2D("wpi0g Theta vs E for photons, only 3g","#theta_{#gamma}","E_{#gamma}",TheBins,EnergyMBins,"O3g_gThevsE",true);
    hO3g_gThevsPhi = hfPi0gH->makeTH2D("wpi0g Theta vs Phi for photons, only 3g","#theta_{#gamma}","#phi_{#gamma}",TheBins,PhiBins,"O3g_gThevsPhi",true);
    hO3g_wThevsE = hfPi0gH->makeTH2D("wpi0g Theta vs E for w, only 3g","#theta_{#omega}","E_{#omega}",TheBins,EnergyMBins,"O3g_wThevsE",true);
    hO3g_wThevsPhi = hfPi0gH->makeTH2D("wpi0g Theta vs Phi for w, only 3g","#theta_{#omega}","#phi_{#omega}",TheBins,PhiBins,"O3g_wThevsPhi",true);
    hO3gOCB_IMgggVsTagCh = hfPi0gH->makeTH2D("wpi0g IMggg vs Tagger channel, only 3g, only CB","IM(#gamma#gamma#gamma)","Tagger channel",IMgggBins,BinSettings(nTagger+20),"O3gOCB_IMgggVsTagCh",true);
    hO3gCBTA_IMgggVsTagCh = hfPi0gH->makeTH2D("wpi0g IMggg vs Tagger channel, only 3g, CB and TAPS","IM(#gamma#gamma#gamma)","Tagger channel",IMgggBins,BinSettings(nTagger+20),"O3gCBTA_IMgggVsTagCh",true);
    hO3gOCB_TagTimeVsChan = hfPi0gH->makeTH2D("wpi0g Tagger time vs channel, only 3g, only CB","Tagger time","Tagger channel",TimeBins1ns,BinSettings(nTagger+20),"hO3gOCB_TagTimeVsChan",true);
    hO3gCBTA_TagTimeVsChan = hfPi0gH->makeTH2D("wpi0g Tagger time vs channel, only 3g, CB and TAPS","Tagger time","Tagger channel",TimeBins1ns,BinSettings(nTagger+20),"hO3gCBTA_TagTimeVsChan",true);
    hO3gOCB_TagCorTimeVsChan = hfPi0gH->makeTH2D("wpi0g Tagger corr. time vs channel, only 3g, only CB","Tagger time - CBsum time","Tagger channel",TimeBins1ns,BinSettings(nTagger+20),"hO3gOCB_TagCorTimeVsChan",true);
    hO3gCBTA_TagCorTimeVsChan = hfPi0gH->makeTH2D("wpi0g Tagger corr. time vs channel, only 3g, CB and TAPS","Tagger time - CBsum time","Tagger channel",TimeBins1ns,BinSettings(nTagger+20),"hO3gCBTA_TagCorTimeVsChan",true);
    hO3gO1p_IMggg = hfPi0gH->makeTH1D("wpi0g IMggg, only 3g 1p","IM(#gamma#gamma#gamma)","",IMgggBins,"wpi0gO3gO1p_IMggg",true);
    hO3gO1p_IMggg_Thggg = hfPi0gH->makeTH2D("wpi0g IMggg vs wTheta, only 3g 1p","IM(#gamma#gamma#gamma","#theta_{#gamma#gamma#gamma}",IMgggBins,TheBins,"wpi0gO3gO1p_IMggg_Thggg",true);
    hO3gO1p_MM = hfPi0gH->makeTH1D("wpi0g MM, only 3g 1p","MM(#gamma#gamma#gamma)","",MMpBins,"wpi0gO3gO1p_MM",true);

}

void scratch_lheijken_gpwppi0g::FillO3gHistos(const str_wpi0gO3g strwpi0g, const TLorentzVector initphoton, const TTaggerHit th, const double ctt, const double tw)
{
    TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());

    if(strwpi0g.onlyinCB){
        hO3gOCB_time->Fill(strwpi0g.avgT3g,tw);
        hO3gOCB_IM->Fill((strwpi0g.LVw).M(),tw);
        hO3gOCB_MM->Fill((initphoton+InitialProtonVec-strwpi0g.LVw).M(),tw);
        hO3gOCB_IMgggVsTagCh->Fill((strwpi0g.LVw).M(),th.Channel,tw);
        hO3gOCB_TagTimeVsChan->Fill(th.Time,th.Channel);
        hO3gOCB_TagCorTimeVsChan->Fill(ctt,th.Channel);
    }
    else {
        hO3gCBTA_time->Fill(strwpi0g.avgT3g,tw);
        hO3gCBTA_IM->Fill((strwpi0g.LVw).M(),tw);
        hO3gCBTA_MM->Fill((initphoton+InitialProtonVec-strwpi0g.LVw).M(),tw);
        hO3gCBTA_IMgggVsTagCh->Fill((strwpi0g.LVw).M(),th.Channel,tw);
        hO3gCBTA_TagTimeVsChan->Fill(th.Time,th.Channel);
        hO3gCBTA_TagCorTimeVsChan->Fill(ctt,th.Channel);
    }
    for(int it=0; it<3; it++){
        hO3g_gThevsE->Fill(strwpi0g.LV3g.at(it).Theta()*radtodeg,strwpi0g.LV3g.at(it).E(),tw);
        hO3g_gThevsPhi->Fill(strwpi0g.LV3g.at(it).Theta()*radtodeg,strwpi0g.LV3g.at(it).Phi()*radtodeg,tw);
    }
    hO3g_wThevsE->Fill(strwpi0g.LVw.Theta()*radtodeg,strwpi0g.LVw.E(),tw);
    hO3g_wThevsPhi->Fill(strwpi0g.LVw.Theta()*radtodeg,strwpi0g.LVw.Phi()*radtodeg,tw);
}

void scratch_lheijken_gpwppi0g::FillRecoCandHistos(const str_recocand rc, const double tw)
{
    for(unsigned int nci=0; nci<rc.NeuCanTheta.size(); nci++){
        hNeuCanThevsCaloE->Fill(rc.NeuCanTheta[nci]*radtodeg,rc.NeuCanCaloE[nci],tw);
        hNeuCanThevsPhi->Fill(rc.NeuCanTheta[nci]*radtodeg,rc.NeuCanPhi[nci]*radtodeg,tw);
        if(rc.NeuCanInCB[nci]){
            hNeuCanCBTime->Fill(rc.NeuCanTime[nci],tw);
            hNeuCanCBCluSize->Fill(rc.NeuCanCluSize[nci],tw);
        }
        else{
            hNeuCanTAPSTime->Fill(rc.NeuCanTime[nci],tw);
            hNeuCanTAPSCluSize->Fill(rc.NeuCanCluSize[nci],tw);
        }
    }
    for(unsigned int cci=0; cci<rc.ChaCanTheta.size(); cci++){
        hChaCanThevsCaloE->Fill(rc.ChaCanTheta[cci]*radtodeg,rc.ChaCanCaloE[cci],tw);
        hChaCanThevsVetoE->Fill(rc.ChaCanTheta[cci]*radtodeg,rc.ChaCanVetoE[cci],tw);
        hChaCanThevsPhi->Fill(rc.ChaCanTheta[cci]*radtodeg,rc.ChaCanPhi[cci]*radtodeg,tw);
        if(rc.ChaCanInCB[cci]){
            hChaCanCBTime->Fill(rc.ChaCanTime[cci],tw);
            hChaCanCBCluSize->Fill(rc.ChaCanCluSize[cci],tw);
        }
        else{
            hChaCanTAPSTime->Fill(rc.ChaCanTime[cci],tw);
            hChaCanTAPSCluSize->Fill(rc.ChaCanCluSize[cci],tw);
        }
    }

}

AUTO_REGISTER_PHYSICS(scratch_lheijken_gpwppi0g)

