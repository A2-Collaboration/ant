#include "Pi0Dalitz.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "utils/Matcher.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/Interpolated.h"
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
    promptrandom(ExpConfig::Setup::Get()),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                             utils::UncertaintyModels::Interpolated::Type_t::MC,
                             make_shared<utils::UncertaintyModels::FitterSergey>())),
    fitter(nullptr, opts->Get<bool>("FitZVertex", true))
{
    fitter.SetZVertexSigma(3.0);

    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    pid_detector = ExpConfig::Setup::GetDetector<expconfig::detector::PID>();
    veto_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPSVeto>();

    CreateHistos();
}


void Pi0Dalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    //-- Check the decay string for MC pattern
    bool isMC=false;
    string decay;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC) && event.MCTrue().ParticleTree){
        decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
        isMC=true;
    }
    else
        decay = "data" + ExpConfig::Setup::Get().GetName();
    vector<bool> WhichMC; // 0 : pi0->eeg, 1 : pi0->gg
    bool MCpi0eeg = false; bool MCpi0gg = false;
    if(isMC){
        MCpi0eeg  = event.MCTrue().ParticleTree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_eeg),
                                                         utils::ParticleTools::MatchByParticleName);
        MCpi0gg  = event.MCTrue().ParticleTree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),
                                                        utils::ParticleTools::MatchByParticleName);
    }
    WhichMC.push_back(MCpi0eeg); WhichMC.push_back(MCpi0gg);

    //-- The reconstructed candidates
    const auto& candidates = event.Reconstructed().Candidates;
    h_AnalysisStat->Fill(1.);

    //-- Set fit model, might later depend on if proton was reconstructed or not
    fitter.SetUncertaintyModel(fit_model);

    //-- Make all possible combinations of the candidates into groups containing one proton and the rest photons
    utils::ProtonPhotonCombs proton_photons(candidates);
    particle_combs_t protphotcombs = proton_photons();

    //-- MC stuff, for the pi0 -> eeg and pi0->gg channels
    vector<TParticlePtr> TruePart;
    vector<TParticlePtr> RecoMatchPart;
    if(isMC){
        for(int i=0; i<nrPartTypes; i++) h_AnalysisStat_RecMat->Fill(0.,(double)i);
        //--- Fetch the true final state particles
        if(MCpi0eeg){
            auto TruePartp_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,event.MCTrue().ParticleTree);
            auto TruePartep_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::ePlus,event.MCTrue().ParticleTree);
            auto TruePartem_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::eMinus,event.MCTrue().ParticleTree);
            auto TruePartg_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Photon,event.MCTrue().ParticleTree);
            TruePart.push_back(TruePartp_eeg); TruePart.push_back(TruePartep_eeg);
            TruePart.push_back(TruePartem_eeg); TruePart.push_back(TruePartg_eeg);
        }
        if(MCpi0gg){
            auto TruePartp_gg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,event.MCTrue().ParticleTree);
            auto TruePartgs_gg = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, event.MCTrue().ParticleTree);
            if(TruePartgs_gg.size() != 2){LOG(INFO)<<"Didn't find two photons in pi0->gg trueMC"; return;}
            TruePart.push_back(TruePartp_gg);
            TruePart.push_back(TruePartgs_gg.front()); TruePart.push_back(TruePartgs_gg.back());
        }

        //--- Match the true particles to the reconstructed candidates
        auto mcparticleslist = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
        const auto mcparticles = mcparticleslist.GetAll();
        DoMatchTrueRecoStuff(mcparticles, TruePart, candidates, RecoMatchPart);

    }

    //-- Loop over the tagger hits
    int cutind;
    for (const auto &tc : event.Reconstructed().TaggerHits) {
        double cortagtime = triggersimu.GetCorrectedTaggerTime(tc);
        promptrandom.SetTaggerTime(cortagtime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        double taggweight = promptrandom.FillWeight();
        TLorentzVector InitialPhotonVec = tc.GetPhotonBeam();

        cutind=0;
        //--- No cuts (except the intrinsic requirement of a prompt taggerhit)
        DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
        DoTriggerStuff(cutind,taggweight);
        DoRecoCandStuff(cutind,candidates,protphotcombs,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
        DoTrueMCStuff(cutind,WhichMC,TruePart,taggweight);

        cutind++;
        //---- CBEsum cut (only on MC, data is always true)
        if(!triggersimu.HasTriggered()) continue;
        DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
        DoTriggerStuff(cutind,taggweight);
        DoRecoCandStuff(cutind,candidates,protphotcombs,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
        DoTrueMCStuff(cutind,WhichMC,TruePart,taggweight);

        cutind++;
        //---- Only allow tagger hits with E < 450 MeV
        if(InitialPhotonVec.E() > 450.) continue;
        DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
        DoTriggerStuff(cutind,taggweight);
        DoRecoCandStuff(cutind,candidates,protphotcombs,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
        DoTrueMCStuff(cutind,WhichMC,TruePart,taggweight);

        cutind++;
        //---- Exactly three reconstructed candidates
        if(candidates.size()==3){
            DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
            DoTriggerStuff(cutind,taggweight);
            DoRecoCandStuff(cutind,candidates,protphotcombs,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
            DoTrueMCStuff(cutind,WhichMC,TruePart,taggweight);

        }

        cutind++;
        //---- At least four reconstructed candidates
        if(candidates.size()>3){
            DoTaggerStuff(cutind,InitialPhotonVec,tc.Time,taggweight);
            DoTriggerStuff(cutind,taggweight);
            DoRecoCandStuff(cutind,candidates,protphotcombs,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
            DoTrueMCStuff(cutind,WhichMC,TruePart,taggweight);
            //----- Do the kinfit on the hypothesis of measured proton (except its energy) and 3 "photons"
            DoKinFitStuff(3,protphotcombs,InitialPhotonVec,fitter,taggweight);
        }
    }


}

void Pi0Dalitz::Finish()
{
    //-- After-treatment of statistics histograms
    double stat = h_RecoTrueMatch->GetBinContent(1);
    double eff, efferr;
    if(stat > 0){
        for(int i=1; i<= h_RecoTrueMatch->GetNbinsX(); i++){
            eff = h_RecoTrueMatch->GetBinContent(i)/stat;
            h_RecoTrueMatch->SetBinContent(i,eff);
            h_RecoTrueMatch->SetBinError(i,0); //until I find a better way of dealing with the errors
        }
    }
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
    auto hfOverview = new HistogramFactory("Overview",HistFac,"");
    auto hfKinFit = new HistogramFactory("KinFit",HistFac,"");

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
    h_IMeegTrue = hfTrueMC->makeTH1D("IMeeg True","IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),"h_IMeegTrue",true);
    h_IMggTrue = hfTrueMC->makeTH1D("IMgg True","IM(#gamma#gamma)","",BinSettings(250,0.,1000.),"h_IMggTrue",true);
    h_RecoTrueMatch = hfTrueMC->makeTH1D("Reco-True matches","Particle","Fraction where a matched reco candidate was found",BinSettings(5,-0.5,4.5),"h_RecoTrueMatch",true);
    h_RecoTrueAngle = hfTrueMC->makeTH2D("Angle between matched reco cand and true mc","Particle","opening angle [deg]",BinSettings(5,-0.5,4.5),BinSettings(180,-0.5,44.5),"h_RecoTrueAngle",true);
    string particlename[] = {"p","eplus","eminus","photon"};
    h_RecoTrueMatch->GetXaxis()->SetBinLabel(1,"# events");
    for(int i=0; i<nrPartTypes; i++){
        h_RecoTrueMatch->GetXaxis()->SetBinLabel(i+2,particlename[i].c_str());
        h_RecoTrueAngle->GetXaxis()->SetBinLabel(i+2,particlename[i].c_str());
        h_EktrueEkrec[i] = hfTrueMC->makeTH1D(Form("Ek_{true}/Ek_{rec} %s",particlename[i].c_str()),"Ek_{true}/Ek_{rec}","",BinSettings(100,0.5,1.5),Form("h_EktrueEkrec%s",particlename[i].c_str()),true);
    }
    h_EktrueEkrec_gg = hfTrueMC->makeTH1D("Ek_{true}/Ek_{rec} for #gamma pairs","Ek_{true}/Ek_{rec}","",BinSettings(100,0.5,1.5),"h_EktrueEkrec_gg",true);
    //--- Reconstructed candidates
    h_PIDMultUsed = hfCandChecks->makeTH2D("Multiple usage of a PID element","Nr times used/event","Element nr",BinSettings(10,-0.5,9.5),BinSettings(npidch),"h_PIDMultUsed",true);
    h_VetoMultUsed = hfCandChecks->makeTH2D("Multiple usage of a Veto element","Nr times used/event","Element nr",BinSettings(10,-0.5,9.5),BinSettings(nvetoch),"h_VetoMultUsed",true);

    string cutname[] = {"NoCuts","CBE","Eg","3Cand","4Cand"};
    string cuttitle[] = {"no cuts","CBEsum","Eg","3 candidates","4 candidates"};
    vector<HistogramFactory> hfTrigChecksCuts, hfTaggChecksCuts, hfCandChecksCuts, hfTrueChecksCuts;
    for(int i=0; i<nrSel; i++){
        //---- Create subfolders for each cut
        auto hfTrigChecksCutTemp = new HistogramFactory(cutname[i],*hfTrigChecks,"");
        hfTrigChecksCuts.push_back(*hfTrigChecksCutTemp);
        auto hfTaggChecksCutTemp = new HistogramFactory(cutname[i],*hfTaggChecks,"");
        hfTaggChecksCuts.push_back(*hfTaggChecksCutTemp);
        auto hfCandChecksCutTemp = new HistogramFactory(cutname[i],*hfCandChecks,"");
        hfCandChecksCuts.push_back(*hfCandChecksCutTemp);
        auto hfMatchedTemp = new HistogramFactory("ForMCMatched",hfCandChecksCuts.at(i),"");
        auto hfTrueChecksCutTemp = new HistogramFactory(cutname[i],*hfTrueMC,"");
        hfTrueChecksCuts.push_back(*hfTrueChecksCutTemp);
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
        h_ThetavsEnergy[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Theta vs energy %s",cuttitle[i].c_str()),"#theta [deg]","Cluster energy",BinSettings(180,-0.5,179.5),BinSettings(200,0.,600),Form("h_ThetavsEnergy%s",cutname[i].c_str()),true);
        for(int j=0; j<nrPartTypes; j++){
            h_EnergyCBvsPID_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("Energy CB vs PID %s for %s",cuttitle[i].c_str(),particlename[j].c_str()),"CB Energy","PID Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyCBvsPID_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_EnergyTAPSvsTVeto_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("Energy TAPS vs TAPSVeto %s for %s",cuttitle[i].c_str(),particlename[j].c_str()),"TAPS Energy","TAPSVeto Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyTAPSvsTAPSVeto_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_ThetavsEnergy_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("Theta vs energy %s for %s",cuttitle[i].c_str(),particlename[j].c_str()),"#theta [deg]","Cluster energy",BinSettings(180,-0.5,179.5),BinSettings(200,0.,600),Form("h_ThetavsEnergy_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_CBEvsNrCr_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("CB E vs NrCr %s for %s",cuttitle[i].c_str(),particlename[j].c_str()),"CB Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_CBEvsNrCr_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_TAPSEvsNrCr_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("TAPS E vs NrCr %s for %s",cuttitle[i].c_str(),particlename[j].c_str()),"TAPS Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_TAPSEvsNrCr_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
        }
        h_EnergyCBvsPID_RecMat[i][nrPartTypes] = hfMatchedTemp->makeTH2D(Form("Energy CB vs PID %s for photon, only 1 PID hit",cuttitle[i].c_str()),"CB Energy","PID Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyCBvsPID_RecMat%sphoton1PIDhit",cutname[i].c_str()),true);
        h_OpAngpphReco_RecMat[i] = hfMatchedTemp->makeTH1D(Form("OpAngpph %s",cuttitle[i].c_str()),"#theta(p - #sum#gamma)","",BinSettings(90,0.,180.),Form("h_OpAngpphReco_RecMat%s",cutname[i].c_str()),true);
        h_NrRecCand[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Nr reconstructed candidates %s",cuttitle[i].c_str()),"Nr cand in CB","Nr cand in TAPS",BinSettings(20,-0.5,19.5),BinSettings(20,-0.5,19.5),Form("h_NrRecCand%s",cutname[i].c_str()),true);
        //---- Histos: True MC
        h_ee_angle[i] = hfTrueChecksCuts.at(i).makeTH1D(Form("Opening angle e^{+}e^{-}, %s",cuttitle[i].c_str()),"Opening angle e^{+}e^{-}","",BinSettings(90,0.,180.),Form("h_ee_angle_%s",cutname[i].c_str()),true);
        for(int j=0; j<nrPartTypes; j++){
            h_ThetavsEnergy_MCTrue[i][j] = hfTrueChecksCuts.at(i).makeTH2D(Form("Energy vs #theta for %s %s",particlename[j].c_str(),cuttitle[i].c_str()),"#theta [deg]","Energy",BinSettings(180,-0.5,179.5),BinSettings(200,0.,600.),Form("h_ThetavsEnergy_MCTrue%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
        }
        //---- delete the temporary histogramfactories
        delete hfTrigChecksCutTemp; delete hfTaggChecksCutTemp; delete hfCandChecksCutTemp; delete hfMatchedTemp; delete hfTrueChecksCutTemp;
    }
    //--- Overview
    h_AnalysisStat = hfOverview->makeTH1D("Analysis statistics","","",BinSettings(3+nrSel,-0.5,(3+nrSel)-0.5),"h_AnalysisStat",true);
    h_AnalysisStat_RecMat = hfOverview->makeTH2D("Analysis statistics","","",BinSettings(3+nrSel,-0.5,(3+nrSel)-0.5),BinSettings(nrPartTypes,-0.5,nrPartTypes-0.5),"h_AnalysisStat_RecMat",true);
    h_AnalysisStat->GetXaxis()->SetBinLabel(2,"Event"); h_AnalysisStat->GetXaxis()->SetBinLabel(3,"Valid tagghit");
    h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(1,"Event"); h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(2,"RecMat");
    h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(3,"Valid tagghit");
    for(int i=1; i<nrSel; i++){
        h_AnalysisStat->GetXaxis()->SetBinLabel(i+3,cutname[i].c_str());
        h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(i+3,cutname[i].c_str());
    }
    for(int i=0; i<nrPartTypes; i++)
        h_AnalysisStat_RecMat->GetYaxis()->SetBinLabel(i+1,particlename[i].c_str());
    //--- KinFit
    //---- 1proton 3photon
    auto hfKinFit_1p3ph = new HistogramFactory("KF1p3ph",*hfKinFit,"");
    h_KF1p3ph_MMcut = hfKinFit_1p3ph->makeTH1D("KF 1p3ph MM(p) cut for p-#gamma combs","MM(p) [MeV]","",BinSettings(250,400.,1400.),"h_KF1p3ph_MMcut",true);
    h_KF1p3ph_EP = hfKinFit_1p3ph->makeTH2D("KF 1p2ph E vs P","E_{initial}^{fit} - E_{final}^{fit}","P_{initial}^{fit} - P_{final}^{fit}",BinSettings(101,-200,200),BinSettings(101,-200,200),"h_KF1p3ph_EP",true);
    for(int i=0; i<8; i++){
        h_KF1p3ph_Prob[i] = hfKinFit_1p3ph->makeTH1D(Form("KF 1p3ph Probability pcut%d",i),"P(#chi^{2})","",BinSettings(1000,0.,1.),Form("h_KF1p3ph_Prob_pcut%d",i),true);
        h_KF1p3ph_Zv[i] = hfKinFit_1p3ph->makeTH1D(Form("KF 1p3ph Z-vertex pcut%d",i),"Z vertex","",BinSettings(100,-15.,15.),Form("h_KF1p3ph_Zv_pcut%d",i),true);
        h_KF1p3ph_IM3g[i] = hfKinFit_1p3ph->makeTH1D(Form("KF 1p3ph IM(e^{+}e^{-}#gamma) pcut%d",i),"IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),Form("h_KF1p3ph_IM3g_pcut%d",i),true);
    }
}

void Pi0Dalitz::DoTrueMCStuff(const int cut, const std::vector<bool> &WhichMC, const std::vector<TParticlePtr> &trueparts, const double &tw)
{
    TParticlePtr pep; TParticlePtr pem;
    for(auto& truepart: trueparts){
        if(truepart->Type() == ParticleTypeDatabase::Proton) h_ThetavsEnergy_MCTrue[cut][en_p]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
        if(truepart->Type() == ParticleTypeDatabase::ePlus){
            h_ThetavsEnergy_MCTrue[cut][en_ep]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
            pep = truepart;
        }
        if(truepart->Type() == ParticleTypeDatabase::eMinus){
            h_ThetavsEnergy_MCTrue[cut][en_em]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
            pem = truepart;
        }
        if(truepart->Type() == ParticleTypeDatabase::Photon) h_ThetavsEnergy_MCTrue[cut][en_g]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
    }
    //-- pi0->eeg
    if(WhichMC.at(0)){
        if(cut == 0) h_IMeegTrue->Fill((*trueparts.at(1)+*trueparts.at(2)+*trueparts.at(3)).M(),tw);
        h_ee_angle[cut]->Fill(TParticle::CalcAngle(pep,pem)*radtodeg,tw);
    }
    //-- pi0->gg
    if(WhichMC.at(1)){
        if(cut == 0) h_IMggTrue->Fill((*trueparts.at(1)+*trueparts.at(2)).M(),tw);
    }


}

void Pi0Dalitz::DoMatchTrueRecoStuff(const TParticleList &allmcpart, const std::vector<TParticlePtr> &trueparts, const TCandidateList &recocands, std::vector<TParticlePtr> &matchrecopart)
{
    const auto matched  = utils::match1to1(allmcpart, recocands.get_ptr_list(),
                                           [] (const TParticlePtr& p1, const TCandidatePtr& p2) {return p1->Angle(*p2);},
                                           {0.0, std_ext::degree_to_radian(15.0)});

    h_RecoTrueMatch->Fill(0);
    int nrgamma = 0; double Ekrecg=0; double Ektrueg=0;
    for(auto& truepart: trueparts){
        TCandidatePtr match = utils::FindMatched(matched,truepart);
        if(match){
            if(truepart->Type() == ParticleTypeDatabase::Proton){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::Proton,match));
                h_RecoTrueMatch->Fill(en_p+1);
                h_RecoTrueAngle->Fill(en_p+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_p]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_p);
            }
            if(truepart->Type() == ParticleTypeDatabase::ePlus){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::ePlus,match));
                h_RecoTrueMatch->Fill(en_ep+1);
                h_RecoTrueAngle->Fill(en_ep+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_ep]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_ep);
            }
            if(truepart->Type() == ParticleTypeDatabase::eMinus){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::eMinus,match));
                h_RecoTrueMatch->Fill(en_em+1);
                h_RecoTrueAngle->Fill(en_em+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_em]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_em);
            }
            if(truepart->Type() == ParticleTypeDatabase::Photon){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,match));
                h_RecoTrueMatch->Fill(en_g+1);
                h_RecoTrueAngle->Fill(en_g+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_g]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_g);
                nrgamma++;
                Ekrecg+=matchrecopart.back()->Ek();
                Ektrueg+=truepart->Ek();
            }
        }
    }
    if(nrgamma==2 && Ekrecg>0)
        h_EktrueEkrec_gg->Fill(Ektrueg/Ekrecg);
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

void Pi0Dalitz::DoRecoCandStuff(const int cut, const TCandidateList &recocands, particle_combs_t ppcomb, const std::vector<TParticlePtr> &recmatparts, const std::vector<bool> &WhichMC, const TLorentzVector &ig, const double &tw)
{
    h_AnalysisStat->Fill(2+cut,tw);

    //-- Fill the reconstructed candidate information
    std::map<int, int> PIDElFreq;
    std::map<int, int> TVetoElFreq;
    int nrcandCB=0; int nrcandTAPS=0;
    for(const auto& currcand : recocands.get_iter()) {
        TClusterPtr caloclu = currcand->FindCaloCluster();
        TClusterPtr vetoclu = nullptr; if(currcand->VetoEnergy) vetoclu = currcand->FindVetoCluster();
        //--- PID/Veto timings, only once per veto hit, i.e. exclude duplicates
        //    also check which pid and tapsveto clusters are used multiple times
        //    two maps are created which keeps track of how often a veto element is used
        if(currcand->VetoEnergy){
            int elnr = vetoclu->CentralElement;
            if(vetoclu->DetectorType == Detector_t::Type_t::PID){
                auto result = PIDElFreq.insert(std::pair<int, int>(elnr, 1));
                if (result.second == false)
                    result.first->second++;
                else
                    h_PIDEvsT[cut]->Fill(vetoclu->Energy,vetoclu->Time,tw);
            }
            else{
                auto result = TVetoElFreq.insert(std::pair<int, int>(elnr, 1));
                if (result.second == false)
                    result.first->second++;
                else
                    h_TVetoEvsT[cut]->Fill(vetoclu->Energy,vetoclu->Time,tw);
            }
        }
        //--- CB related info
        if(caloclu->DetectorType == Detector_t::Type_t::CB){
            nrcandCB++;
            h_CBEvsT[cut]->Fill(currcand->CaloEnergy,caloclu->Time,tw);
            if(currcand->VetoEnergy && vetoclu->DetectorType == Detector_t::Type_t::PID){
                h_EnergyCBvsPID[cut]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
                h_TimeCBvsPID[cut]->Fill(currcand->Time,vetoclu->Time,tw);
                h_ThPhCBvsPID[cut]->Fill((caloclu->Position.Theta()-vetoclu->Position.Theta())*radtodeg,(caloclu->Position.Phi()-vetoclu->Position.Phi())*radtodeg,tw);
                h_CBEvsNrCr[cut]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
            }
        }
        //--- TAPS related info
        else{
            nrcandTAPS++;
            h_TAPSEvsT[cut]->Fill(currcand->CaloEnergy,currcand->Time,tw);
            //---- Here I'm ignoring the possibility (is it even possible, I don't think so) that a TAPS cluster can be paired with a PID cluster
            if(currcand->VetoEnergy && vetoclu->DetectorType == Detector_t::Type_t::TAPSVeto){
                h_EnergyTAPSvsTVeto[cut]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
                h_TimeTAPSvsTVeto[cut]->Fill(currcand->Time,vetoclu->Time,tw);
                h_ThPhTAPSvsTVeto[cut]->Fill((caloclu->Position.Theta()-vetoclu->Position.Theta())*radtodeg,(caloclu->Position.Phi()-vetoclu->Position.Phi())*radtodeg,tw);
                h_TAPSEvsNrCr[cut]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
            }
        }
        h_ThetavsEnergy[cut]->Fill(currcand->Theta*radtodeg,currcand->CaloEnergy,tw);
    }
    h_NrRecCand[cut]->Fill(nrcandCB,nrcandTAPS,tw);

    //-- How often are each PID/Veto element used per event
    if(cut==0){
        for(auto elem : PIDElFreq)
            h_PIDMultUsed->Fill(elem.second,elem.first,tw);
        for(auto elem : TVetoElFreq)
            h_VetoMultUsed->Fill(elem.second,elem.first,tw);
    }

    //-- Fill the info from the proton-photon particle combinations
    TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
    TLorentzVector InitialStateVec = InitialProtonVec + ig;
    TVector3 InStBoost = -InitialStateVec.BoostVector();
    if(ppcomb.size()>0){
        //--- Fill combinatorially the MM(p) and opening angle between the proton and the photons in CM
        //    for this, at least one proton and one photon is needed
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

    //-- Fill histograms for the candidates mathed with true MC
    TLorentzVector matvec_p;
    vector<TLorentzVector> matvec_pi0;
    int pt;
    for(auto& recmatpart: recmatparts){
        //--- What kind of particle is it
        if(recmatpart->Type() == ParticleTypeDatabase::Proton){pt = en_p; matvec_p = *recmatpart;}
        else if(recmatpart->Type() == ParticleTypeDatabase::ePlus){ pt = en_ep; matvec_pi0.push_back(*recmatpart);}
        else if(recmatpart->Type() == ParticleTypeDatabase::eMinus){ pt = en_em; matvec_pi0.push_back(*recmatpart);}
        else if (recmatpart->Type() == ParticleTypeDatabase::Photon){ pt = en_g; matvec_pi0.push_back(*recmatpart);}
        else continue;
        h_AnalysisStat_RecMat->Fill(2+cut,pt,tw);
        TCandidatePtr currcand = recmatpart->Candidate;
        h_ThetavsEnergy_RecMat[cut][pt]->Fill(radtodeg*currcand->Theta,currcand->CaloEnergy,tw);
        //--- CB related info
        if(currcand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_CBEvsNrCr_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
            if(currcand->VetoEnergy && currcand->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_EnergyCBvsPID_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
                if(PIDElFreq.at(currcand->FindVetoCluster()->CentralElement) == 1 && pt == en_g)
                    h_EnergyCBvsPID_RecMat[cut][nrPartTypes]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
            }
        }
        //--- TAPS related info
        else{
            h_TAPSEvsNrCr_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
            if(currcand->VetoEnergy && currcand->FindVetoCluster()->DetectorType == Detector_t::Type_t::TAPSVeto){
                h_EnergyTAPSvsTVeto_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
            }
        }
    }
    //--- Boost to CM frame and check the proton-pi0 opening angle
    if(((matvec_pi0.size()==3 && WhichMC.at(0)) || (matvec_pi0.size()==2 && WhichMC.at(1))) && matvec_p.E()){
        TLorentzVector summatvec_pi0;
        for(auto& matvec: matvec_pi0) summatvec_pi0 += matvec;
        summatvec_pi0.Boost(InStBoost);
        matvec_p.Boost(InStBoost);
        h_OpAngpphReco_RecMat[cut]->Fill(matvec_p.Angle(summatvec_pi0.Vect())*radtodeg,tw);
    }
}

void Pi0Dalitz::DoKinFitStuff(const int nrph, particle_combs_t ppcomb, const TLorentzVector &ig, utils::KinFitter &fitobj, const double tw)
{
    TParticleList bestfitrecph, bestfitph;
    TParticlePtr bestfitrecp, bestfitp;
    double fit_z_vert = -50;
    double fitbeamE = -50;
    double fitprob = -50;
    //-- loop over the proton-photon combinations
    for(const auto& comb : ppcomb) {
        TParticlePtr protontofit = comb.Proton;
        //-- select on the number of photons required in the current kinfit hypothesis
        for(auto photoncomb = analysis::utils::makeCombination(comb.Photons,nrph); !photoncomb.done(); ++photoncomb ){
            TParticleList photonstofit;
            TLorentzVector photonssum;
            for(int i=0; i<nrph; i++){
                photonstofit.push_back(photoncomb.at(i));
                photonssum += *photoncomb.at(i);
            }
            //-- check that the MM(p) is near mp
            double mmp = (LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass()) + ig - photonssum).M();
            if(!ParticleTypeDatabase::Proton.GetWindow(350).Round().Contains(mmp))
                continue;
            h_KF1p3ph_MMcut->Fill(mmp,tw);
            //-- do the kinfit
            const auto& result = fitobj.DoFit(ig.E(), protontofit, photonstofit);
            //-- if the fit did not converge or did not result in a better probability, try next combination
            if(result.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(fitprob, result.Probability))
                continue;

            //-- else, store the fit result and the reconstructed input
            bestfitrecp = protontofit;
            bestfitrecph = photonstofit;
            bestfitp = fitobj.GetFittedProton();
            bestfitph = fitobj.GetFittedPhotons();
            fit_z_vert = fitobj.GetFittedZVertex();
            fitbeamE = fitobj.GetFittedBeamE();
        }
    }

    if(fitprob==(-50)) return;
    //-- fill histograms
    TLorentzVector fitphsum; for(int i=0; i<nrph; i++) fitphsum = fitphsum + *bestfitph.at(i);
    TLorentzVector fitinitial = LorentzVec({0, 0, fitbeamE}, fitbeamE) + LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass());
    TLorentzVector IFdiff = fitinitial - *bestfitp - fitphsum;
    h_KF1p3ph_EP->Fill(IFdiff.E(),IFdiff.P());
    double pcuts[] = {0., 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.5};
    for(int i=0; i<8; i++){
        if(fitprob<pcuts[i]) break;
        h_KF1p3ph_Prob[i]->Fill(fitprob,tw);
        h_KF1p3ph_Zv[i]->Fill(fit_z_vert,tw);
        h_KF1p3ph_IM3g[i]->Fill(fitphsum.M(),tw);
    }
}

AUTO_REGISTER_PHYSICS(Pi0Dalitz)

