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
    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    pid_detector = ExpConfig::Setup::GetDetector<expconfig::detector::PID>();
    veto_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPSVeto>();

    CreateHistos();
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

    //-- Collect reconstructed info for each candidate
    const auto& candidates = event.Reconstructed().Candidates;
    //--- Store energy, angles, time, and other info in vectors depending on their data type
    //    this will be collected in a reco info struct
    vector<RecoCandInfo> RecoCandInfos;
    vector<double> doublelist;
    vector<int> intlist;
    vector<bool> boollist;
    std::map<int, int> PIDElFreq;
    std::map<int, int> TVetoElFreq;
    for(const auto& cand : candidates.get_iter()) {
        doublelist.push_back(cand->CaloEnergy);
        doublelist.push_back(cand->VetoEnergy);
        doublelist.push_back(cand->Theta);
        doublelist.push_back(cand->Phi);
        doublelist.push_back(cand->Time);
        intlist.push_back(cand->ClusterSize);
        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
            boollist.push_back(true);
        else
            boollist.push_back(false);
        auto vetoclu = cand->FindVetoCluster();
        if(vetoclu){
            doublelist.push_back((vetoclu->Position).Theta());
            doublelist.push_back((vetoclu->Position).Phi());
            doublelist.push_back(vetoclu->Time);
            intlist.push_back(vetoclu->CentralElement);
            if(vetoclu->DetectorType == Detector_t::Type_t::PID)
                boollist.push_back(true);
            else
                boollist.push_back(false);
        }
        else {
            doublelist.push_back(-10000);
            doublelist.push_back(-10000);
            doublelist.push_back(-10000);
            intlist.push_back(-10000);
            boollist.push_back(false);
        }
        //--- Check which pid and tapsveto clusters are used multiple times
        //    VetoElClean contains only the veto el nr if it is the first occurance
        int vetoelclean = intlist.at(1);
        CheckVetoUsage(intlist.at(1), boollist.at(1), vetoelclean, PIDElFreq, TVetoElFreq);
        intlist.push_back(vetoelclean);
        //--- Create and store the RecoInfo struct for this candidate
        RecoCandInfo rcand;
        FillRecoInfo(&rcand, doublelist, intlist, boollist);
        RecoCandInfos.push_back(rcand);
        doublelist.clear(); intlist.clear(); boollist.clear();
    }

    //-- Make all possible combinations of the candidates into groups containing one proton and the rest photons
    utils::ProtonPhotonCombs proton_photons(candidates);
    particle_combs_t protphotcombs = proton_photons();

    //-- Loop over the tagger hits
    int cutind;
    for (const auto &tc : event.Reconstructed().TaggerHits) {
        double cortagtime = triggersimu.GetCorrectedTaggerTime(tc);
        promptrandom.SetTaggerTime(cortagtime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        double taggweight = promptrandom.FillWeight();
        TLorentzVector InitialPhotonVec = tc.GetPhotonBeam();

        //--- Do some stuff with the true MC info
        if(MCpi0eeg)
            DoTrueMCStuff(1, TrueVecep, TrueVecem, TrueVecGammas);

        cutind=0;
        //--- No cuts (except the intrinsic requirement of a prompt taggerhit)
        DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
        DoTriggerStuff(cutind,taggweight);
        DoRecoCandStuff(cutind,RecoCandInfos,PIDElFreq,TVetoElFreq,protphotcombs,InitialPhotonVec,taggweight);

        cutind++;
        //---- CBEsum cut (only on MC, data is always true)
        if(!triggersimu.HasTriggered()) continue;
        DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
        DoTriggerStuff(cutind,taggweight);
        DoRecoCandStuff(cutind,RecoCandInfos,PIDElFreq,TVetoElFreq,protphotcombs,InitialPhotonVec,taggweight);

        cutind++;
        //---- Only allow tagger hits with E < 450 MeV
        if(InitialPhotonVec.E() > 450.) continue;
        DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
        DoTriggerStuff(cutind,taggweight);
        DoRecoCandStuff(cutind,RecoCandInfos,PIDElFreq,TVetoElFreq,protphotcombs,InitialPhotonVec,taggweight);

        cutind++;
        //---- Exactly three reconstructed candidates
        if(RecoCandInfos.size()==3){
            DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
            DoTriggerStuff(cutind,taggweight);
            DoRecoCandStuff(cutind,RecoCandInfos,PIDElFreq,TVetoElFreq,protphotcombs,InitialPhotonVec,taggweight);
        }

        cutind++;
        //---- Exactly four reconstructed candidates
        if(RecoCandInfos.size()==4){
            DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
            DoTriggerStuff(cutind,taggweight);
            DoRecoCandStuff(cutind,RecoCandInfos,PIDElFreq,TVetoElFreq,protphotcombs,InitialPhotonVec,taggweight);
        }
    }


}

void Pi0Dalitz::Finish()
{
}

void Pi0Dalitz::ShowResult()
{

}

void Pi0Dalitz::CreateHistos()
{
    //-- Histogram folders
    auto hfTrueMC = new HistogramFactory("TrueMC",HistFac,"");
    auto hfTrigChecks = new HistogramFactory("TrigChecks",HistFac,"");
    auto hfTaggChecks = new HistogramFactory("TaggChecks",HistFac,"");
    auto hfCandChecks = new HistogramFactory("CandChecks",HistFac,"");

    //-- Specific bin settings
    //--- Tagger: make some bin edges half way between the energy values.
    vector<double> photonenergies;
    int ntaggch = tagger_detector->GetNChannels();
    for(int ch=0;ch<ntaggch;ch++)
        photonenergies.push_back(tagger_detector->GetPhotonEnergy(ch));
    vector<double> energy_bins(ntaggch+1);
    for(int b=0;b<(ntaggch-1);b++){
      energy_bins.at(b+1)=0.5*(photonenergies.at(b)+photonenergies.at(b+1));
    }
    //--- Tagger: and the top and bottom have width of the adjacent bin
    energy_bins.at(0) = energy_bins.at(1) - (energy_bins.at(2) - energy_bins.at(1));
    energy_bins.at(ntaggch) = energy_bins.at(ntaggch-1) + (energy_bins.at(ntaggch-1)- energy_bins.at(ntaggch-2));

    int npidch = pid_detector->GetNChannels();
    int nvetoch = veto_detector->GetNChannels();

    //-- The histograms
    //--- MC true
    h_IsPi0eegMC = hfTrueMC->makeTH1D("Is Pi0eeg","","",BinSettings(10,-0.5,9.5),"h_IsPi0eeg",true);
    h_IMeegTrue = hfTrueMC->makeTH1D("IMeeg True","IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),"h_IMeegTrue",true);
    //--- Reconstructed candidates
    h_PIDMultUsed = hfCandChecks->makeTH2D("Multiple usage of a PID element","Nr times used/event","Element nr",BinSettings(10,-0.5,9.5),BinSettings(npidch),"h_PIDMultUsed",true);
    h_VetoMultUsed = hfCandChecks->makeTH2D("Multiple usage of a Veto element","Nr times used/event","Element nr",BinSettings(10,-0.5,9.5),BinSettings(nvetoch),"h_VetoMultUsed",true);

    string cutname[] = {"NoCuts","CBE","Eg","3Cand","4Cand"};
    string cuttitle[] = {"no cuts","CBEsum","Eg","3 candidates","4 candidates"};
    vector<HistogramFactory> hfTrigChecksCuts, hfTaggChecksCuts, hfCandChecksCuts;
    for(int i=0; i<nrSel; i++){
        //---- Create subfolders for each cut
        auto hfTrigChecksCutTemp = new HistogramFactory(cutname[i],*hfTrigChecks,"");
        hfTrigChecksCuts.push_back(*hfTrigChecksCutTemp);
        auto hfTaggChecksCutTemp = new HistogramFactory(cutname[i],*hfTaggChecks,"");
        hfTaggChecksCuts.push_back(*hfTaggChecksCutTemp);
        auto hfCandChecksCutTemp = new HistogramFactory(cutname[i],*hfCandChecks,"");
        hfCandChecksCuts.push_back(*hfCandChecksCutTemp);
        //---- Histos: Trigger info
        h_CBEsum[i] = hfTrigChecksCuts.at(i).makeTH1D(Form("CBEsum %s",cuttitle[i].c_str()),"CB Esum","",BinSettings(250,0.,2000.),Form("h_CBEsum%s",cutname[i].c_str()),true);
        //---- Histos: Tagger info
        h_TaggTimeuw[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger time  unweighted%s",cuttitle[i].c_str()),"Time","",BinSettings(200,-400,400),Form("h_TaggTimeuw%s",cutname[i].c_str()),true);
        h_TaggTimeww[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger time  with weights%s",cuttitle[i].c_str()),"Time","",BinSettings(200,-400,400),Form("h_TaggTimeww%s",cutname[i].c_str()),true);
        h_TaggPhEnww[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger photon energy with weights %s",cuttitle[i].c_str()),energy_bins,"Photon energy","",Form("h_TaggPhEnww%s",cutname[i].c_str()),true);
        h_TaggPhEnuw[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger photon energy unweighted %s",cuttitle[i].c_str()),energy_bins,"Photon energy","",Form("h_TaggPhEnuw%s",cutname[i].c_str()),true);
        //---- Histos: Reconstructed candidates
        h_IMeegReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("IMeeg %s",cuttitle[i].c_str()),"IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),Form("h_IMeegReco%s",cutname[i].c_str()),true);
        h_IMggReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("IMgg %s",cuttitle[i].c_str()),"IM(#gamma#gamma)","",BinSettings(250,0.,1000.),Form("h_IMggReco%s",cutname[i].c_str()),true);
        h_MMpReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("MMp %s",cuttitle[i].c_str()),"MM(p)","",BinSettings(250,400.,1400.),Form("h_MMpReco%s",cutname[i].c_str()),true);
        h_OpAngpphReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("OpAngpph %s",cuttitle[i].c_str()),"#theta(p - #sum#gamma)","",BinSettings(90,0.,180.),Form("h_OpAngpphReco%s",cutname[i].c_str()),true);
        h_PIDEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("PID E vs T %s",cuttitle[i].c_str()),"PID Energy","PID Time",BinSettings(200,0.,20),BinSettings(241,-60.5,60.5),Form("h_PIDEvsT%s",cutname[i].c_str()),true);
        h_TVetoEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("TAPSVeto E vs T %s",cuttitle[i].c_str()),"TAPSVeto Energy","TAPSVeto Time",BinSettings(200,0.,20),BinSettings(241,-60.5,60.5),Form("h_TVetoEvsT%s",cutname[i].c_str()),true);
        h_CBEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("CB E vs T %s",cuttitle[i].c_str()),"CB Energy","CB Time",BinSettings(200,0.,600),BinSettings(241,-60.5,60.5),Form("h_CBEvsT%s",cutname[i].c_str()),true);
        h_TAPSEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("TAPS E vs T %s",cuttitle[i].c_str()),"TAPS Energy","TAPS Time",BinSettings(200,0.,600),BinSettings(241,-60.5,60.5),Form("h_TAPSEvsT%s",cutname[i].c_str()),true);
        h_CBEvsNrCr[i] = hfCandChecksCuts.at(i).makeTH2D(Form("CB E vs NrCr %s",cuttitle[i].c_str()),"CB Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_CBEvsNrCr%s",cutname[i].c_str()),true);
        h_TAPSEvsNrCr[i] = hfCandChecksCuts.at(i).makeTH2D(Form("TAPS E vs NrCr %s",cuttitle[i].c_str()),"TAPS Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_TAPSEvsNrCr%s",cutname[i].c_str()),true);
        h_TimeCBvsPID[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Time CB vs PID %s",cuttitle[i].c_str()),"CB Time","PID Time",BinSettings(241,-60.5,60.5),BinSettings(241,-60.5,60.5),Form("h_TimeCBvsPID%s",cutname[i].c_str()),true);
        h_EnergyCBvsPID[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Energy CB vs PID %s",cuttitle[i].c_str()),"CB Energy","PID Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyCBvsPID%s",cutname[i].c_str()),true);
        h_TimeTAPSvsTVeto[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Time TAPS vs TAPSVeto %s",cuttitle[i].c_str()),"TAPS Time","TAPSVeto Time",BinSettings(241,-60.5,60.5),BinSettings(241,-60.5,60.5),Form("h_TimeTAPSvsTAPSVeto%s",cutname[i].c_str()),true);
        h_EnergyTAPSvsTVeto[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Energy TAPS vs TAPSVeto %s",cuttitle[i].c_str()),"TAPS Energy","TAPSVeto Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyTAPSvsTAPSVeto%s",cutname[i].c_str()),true);
        h_ThPhCBvsPID[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Theta-Phi CB vs PID %s",cuttitle[i].c_str()),"#theta_{CB-PID}","#phi_{CB-PID}",BinSettings(400,-100.,100.),BinSettings(400,-100.,100.),Form("h_ThPhCBvsPID%s",cutname[i].c_str()),true);
        h_ThPhTAPSvsTVeto[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Theta-Phi angle TAPS vs TAPSVeto %s",cuttitle[i].c_str()),"#theta_{TAPS-TAPSVeto}","#phi_{TAPS-TAPSVeto}",BinSettings(400,-100.,100.),BinSettings(400,-100.,100.),Form("h_ThPhTAPSvsTVeto%s",cutname[i].c_str()),true);
        //---- delete the temporary histogramfactories
        delete hfTrigChecksCutTemp; delete hfTaggChecksCutTemp; delete hfCandChecksCutTemp;
    }
}

void Pi0Dalitz::CheckVetoUsage(const int elnr, const bool inpid, int &elnrcleaned, std::map<int,int> &pidfreq, std::map<int,int> &vetofreq)
{
    //--- Check which pid and tapsveto clusters are used multiple times
    //    two maps are created which keeps track of how often a veto element is used
    //    and a cleaned element vector where only the index of the first occurance of a multiply used veto element is kept
    if(elnr>-1){
        if(inpid){
            auto result = pidfreq.insert(std::pair<int, int>(elnr, 1));
            if (result.second == false){
                result.first->second++;
                elnrcleaned = -10000;
            }
        }
        else{
            auto result = vetofreq.insert(std::pair<int, int>(elnr, 1));
            if (result.second == false){
                result.first->second++;
                elnrcleaned = -10000;
            }
        }
    }

}

void Pi0Dalitz::FillRecoInfo(RecoCandInfo *rc, const std::vector<double> doublelist, const std::vector<int> intlist, const std::vector<bool> boollist)
{
    //  THIS IS STUPID, SHOULD USE MAP INSTEAD

    //-- deposited energies
    rc->CaloE=doublelist.at(0); rc->VetoE=doublelist.at(1);
    //-- reconstructed angles and time for the cluster
    rc->Theta=doublelist.at(2); rc->Phi=doublelist.at(3); rc->CaloT=doublelist.at(4);
    //-- number of crystals in calo cluster
    rc->CaloCluSize=intlist.at(0);
    //-- bool if calo=CB
    rc->InCB=boollist.at(0);

    //-- angles and time for the PID/Veto element
    rc->VetoTheta=doublelist.at(5); rc->VetoPhi=doublelist.at(6); rc->VetoT=doublelist.at(7);
    //-- element number of PID/veto (-10000 if not existing)
    //   and element number of veto if it's the first candidate with this PID/veto assigned (-10000 if not)
    rc->VetoEl=intlist.at(1); rc->VetoElClean=intlist.at(2);
    //-- bool  if veto=PID
    rc->InPID=boollist.at(1);
}


void Pi0Dalitz::DoTrueMCStuff(const int WhichMC, const TLorentzVector &ep, const TLorentzVector &em, const vector<TLorentzVector> &g)
{
    h_IsPi0eegMC->Fill(WhichMC);
    if(g.size()>0)
        h_IMeegTrue->Fill((ep+em+g.at(0)).M());
}

void Pi0Dalitz::DoTaggerStuff(const int cut, const TLorentzVector &g, const double &time, const double &tw)
{
    h_TaggTimeuw[cut]->Fill(time);
    h_TaggTimeww[cut]->Fill(time,tw);
    h_TaggPhEnww[cut]->Fill(g.E(),tw);
    h_TaggPhEnuw[cut]->Fill(g.E());

}

void Pi0Dalitz::DoTriggerStuff(const int cut, const double &tw)
{
    h_CBEsum[cut]->Fill(triggersimu.GetCBEnergySum(),tw);
}

void Pi0Dalitz::DoRecoCandStuff(const int cut, const vector<RecoCandInfo> &recocandinfo, std::map<int,int> &pidfreq, std::map<int,int> &vetofreq, particle_combs_t ppcomb, const TLorentzVector &ig, const double &tw)
{
    //-- How often are each PID/Veto element used per event
    for(auto elem : pidfreq)
        h_PIDMultUsed->Fill(elem.second,elem.first,tw);
    for(auto elem : vetofreq)
        h_VetoMultUsed->Fill(elem.second,elem.first,tw);

    //-- Fill the reconstructed candidate information
    for(unsigned i=0; i<recocandinfo.size();i++){
        //--- PID/Veto timings, only once per veto hit, i.e. exclude duplicates
        if(recocandinfo.at(i).VetoElClean>-1){
            if(recocandinfo.at(i).InPID)
                h_PIDEvsT[cut]->Fill(recocandinfo.at(i).VetoE,recocandinfo.at(i).VetoT,tw);
            else
                h_TVetoEvsT[cut]->Fill(recocandinfo.at(i).VetoE,recocandinfo.at(i).VetoT,tw);
        }
        //--- CB related info
        if(recocandinfo.at(i).InCB){
            h_CBEvsT[cut]->Fill(recocandinfo.at(i).CaloE,recocandinfo.at(i).CaloT,tw);
            if(recocandinfo.at(i).VetoEl>-1 && recocandinfo.at(i).InPID){
                h_EnergyCBvsPID[cut]->Fill(recocandinfo.at(i).CaloE,recocandinfo.at(i).VetoE,tw);
                h_TimeCBvsPID[cut]->Fill(recocandinfo.at(i).CaloT,recocandinfo.at(i).VetoT,tw);
                h_ThPhCBvsPID[cut]->Fill((recocandinfo.at(i).Theta-recocandinfo.at(i).VetoTheta)*radtodeg,(recocandinfo.at(i).Phi-recocandinfo.at(i).VetoPhi)*radtodeg,tw);
                h_CBEvsNrCr[cut]->Fill(recocandinfo.at(i).CaloE,recocandinfo.at(i).CaloCluSize,tw);
            }
        }
        //--- TAPS related info
        else{
            h_TAPSEvsT[cut]->Fill(recocandinfo.at(i).CaloE,recocandinfo.at(i).CaloT,tw);
            //---- Here I'm ignoring the possibility (is it even possible) that a TAPS cluster can be paired with a PID cluster
            if(recocandinfo.at(i).VetoEl>-1 && !(recocandinfo.at(i).InPID)){
                h_EnergyTAPSvsTVeto[cut]->Fill(recocandinfo.at(i).CaloE,recocandinfo.at(i).VetoE,tw);
                h_TimeTAPSvsTVeto[cut]->Fill(recocandinfo.at(i).CaloT,recocandinfo.at(i).VetoT,tw);
                h_ThPhTAPSvsTVeto[cut]->Fill((recocandinfo.at(i).Theta-recocandinfo.at(i).VetoTheta)*radtodeg,(recocandinfo.at(i).Phi-recocandinfo.at(i).VetoPhi)*radtodeg,tw);
                h_TAPSEvsNrCr[cut]->Fill(recocandinfo.at(i).CaloE,recocandinfo.at(i).CaloCluSize,tw);
            }
        }
    }
    //-- Fill the info from the proton-photon particle combinations
    if(ppcomb.size()>0){
        //--- Fill combinatorially the MM(p) and opening angle between the proton and the photons in CM
        //    for this, at least one proton and one photon is needed
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
        TLorentzVector InitialStateVec = InitialProtonVec + ig;
        TVector3 InStBoost = -InitialStateVec.BoostVector();
        for(const auto& comb : ppcomb){
            TLorentzVector TLproton = *comb.Proton;
            TLorentzVector TLallg;
            if(comb.Photons.size() < 1) break;
            for(const auto& ph : comb.Photons)
                TLallg += *ph;
            h_MMpReco[cut]->Fill((InitialProtonVec + ig - TLallg).M(),tw);
            //---- Boost the proton and the summed photons vector to CM and check their opening angle
            TLproton.Boost(InStBoost);
            TLallg.Boost(InStBoost);
            h_OpAngpphReco[cut]->Fill(TLproton.Angle(TLallg.Vect())*radtodeg,tw);
        }
        //-- Assume all candidates are photons and combinatorially fill IM histograms
        auto firstcomb = ppcomb.begin();
        TParticleList AllPhotons = firstcomb->Photons;
        AllPhotons.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, firstcomb->Proton->Candidate));
        //--- First IM(gg)
        if(AllPhotons.size() > 1){
            for(auto photoncomb = analysis::utils::makeCombination(AllPhotons,2); !photoncomb.done(); ++photoncomb ){
                const auto& ph1 = photoncomb.at(0);
                const auto& ph2 = photoncomb.at(1);
                h_IMggReco[cut]->Fill((*ph1 + *ph2).M());
            }
        }
        //--- Then IM(ggg)
        if(AllPhotons.size() > 2){
            for(auto photoncomb = analysis::utils::makeCombination(AllPhotons,3); !photoncomb.done(); ++photoncomb ){
                const auto& ph1 = photoncomb.at(0);
                const auto& ph2 = photoncomb.at(1);
                const auto& ph3 = photoncomb.at(2);
                h_IMeegReco[cut]->Fill((*ph1 + *ph2 + *ph3).M());
            }
        }
    }
}

AUTO_REGISTER_PHYSICS(Pi0Dalitz)

