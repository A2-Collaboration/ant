#include "checkcb.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/Trigger.h"
#include "base/Logger.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "base/WrapTFile.h"

#include <iostream>
#include <memory>

/**
 * checking info from cb
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_lheijken_checkcb::scratch_lheijken_checkcb(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()),
    tdc_converter(expconfig::detector::Trigger::Reference_CATCH_CBCrate)
{
    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    pid_detector = ExpConfig::Setup::GetDetector<expconfig::detector::PID>();

    CreateHistos();

}


void scratch_lheijken_checkcb::ProcessEvent(const TEvent& event, manager_t&)
{

    triggersimu.ProcessEvent(event);

    const double TriggerRefTime = triggersimu.GetRefTiming();
    hTrigRefTiming->Fill(TriggerRefTime);

    //-- Loop over DetectorReadhits in event and extract the ones from CB to look at
    double CBesum = 0;
    int elemhastiming[cb_detector->GetNChannels()];
    memset(elemhastiming, 0, cb_detector->GetNChannels()*sizeof(int));
    const auto& readhits = event.Reconstructed().DetectorReadHits;
    for(const TDetectorReadHit& readhit : readhits) {
        if(readhit.DetectorType != Detector_t::Type_t::CB)
            continue;
        if(readhit.ChannelType == Channel_t::Type_t::Timing) {
            int tmult=0;
            for(auto& time : readhit.Values){
                hDRHUncalTimeAll->Fill(time.Uncalibrated,readhit.Channel);
                hDRHCalTimeAll->Fill(time.Calibrated,readhit.Channel);
                if(tmult==0){
                    hDRHUncalTimeFirst->Fill(time.Uncalibrated,readhit.Channel);
                    hDRHCalTimeFirst->Fill(time.Calibrated,readhit.Channel);
                }
                tmult++;
            }
            hDRHUncalTimeMult->Fill(tmult,readhit.Channel);
            elemhastiming[readhit.Channel] = elemhastiming[readhit.Channel]+1;
        }
        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            int emult=0;
            for(auto& integral : readhit.Values){
                hDRHUncalEnergy->Fill(integral.Uncalibrated,readhit.Channel);
                hDRHCalEnergy->Fill(integral.Calibrated,readhit.Channel);
                emult++;
                CBesum += integral.Calibrated;
            }
            hDRHUncalEnergyMult->Fill(emult,readhit.Channel);
        }
    }
    hDRHCBesum->Fill(CBesum);
    for(const TDetectorReadHit& readhit : readhits) {
        if(readhit.DetectorType != Detector_t::Type_t::CB)
            continue;
        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            for(auto& integral : readhit.Values){
                if(elemhastiming[readhit.Channel]==0){
                    hDRHUnCalEn_woTiming->Fill(integral.Uncalibrated,readhit.Channel);
                    hDRHCalEn_woTiming->Fill(integral.Calibrated,readhit.Channel);
                }
                else {
                    hDRHUnCalEn_wTiming->Fill(integral.Uncalibrated,readhit.Channel);
                    hDRHCalEn_wTiming->Fill(integral.Calibrated,readhit.Channel);
                }
            }
        }
    }
    //-- Loop over Cluster and Clusterhits in event and extract the ones from CB to look at
    int nrCBclu = 0;
    auto& clusters = event.Reconstructed().Clusters;
    for (const TCluster& cluster : clusters){
        if(cluster.DetectorType != Detector_t::Type_t::CB )
            continue;
        for(const TClusterHit& hit : cluster.Hits){
            hCHTime->Fill(hit.Time,hit.Channel);
            hCHEnergy->Fill(hit.Energy,hit.Channel);
//            for(auto& datum : hit.Data) {
//                if(datum.Type == Channel_t::Type_t::Integral) {
//                    hCHTimeRawE->Fill(hit.Time,std::log10(datum.Value.Uncalibrated),hit.Channel);
//                }
//            }            
        }
        hClTime->Fill(cluster.Time,cluster.CentralElement);
        hClEnergy->Fill(cluster.Energy,cluster.CentralElement);
        hClTimeVsEnergy->Fill(cluster.Time,cluster.Energy);
        nrCBclu++;
    }
    //-- Checking PID CB Phi matching
    TClusterPtr cluster_pid;
    TClusterPtr cluster_cb;
    bool onePIDCl = true;
    bool oneCBCl = true;
    for(auto it_cl = clusters.begin(); it_cl != clusters.end(); ++it_cl) {
        if(it_cl->DetectorType == Detector_t::Type_t::PID) {
            if(cluster_pid)
                onePIDCl = false;
            cluster_pid = it_cl.get_ptr();
        }
        else if(it_cl->DetectorType == Detector_t::Type_t::CB) {
            if(cluster_cb)
                oneCBCl = false;
            cluster_cb = it_cl.get_ptr();
        }
    }
    if(cluster_pid && cluster_cb && onePIDCl && oneCBCl) {
        const double phi_cb_degrees = std_ext::radian_to_degree(cluster_cb->Position.Phi());
        const double phi_pid_degrees = std_ext::radian_to_degree(cluster_pid->Position.Phi());
        hCBPhiPIDPhi->Fill(phi_cb_degrees,phi_pid_degrees);
    }
    //-- Also fill the CB esum depending on the number of clusters in CB
    if(nrCBclu==1) hDRHCBesum_1Cl->Fill(CBesum);
    if(nrCBclu==2) hDRHCBesum_2Cl->Fill(CBesum);
    if(nrCBclu==3) hDRHCBesum_3Cl->Fill(CBesum);
    if(nrCBclu==4) hDRHCBesum_4Cl->Fill(CBesum);
    if(nrCBclu==5) hDRHCBesum_5Cl->Fill(CBesum);

    //-- Loop over candidates and look at their CB clusters
    const auto& cands = event.Reconstructed().Candidates;
    for (const TCandidate& cand : cands){
            const auto calclu = cand.FindCaloCluster();
            if(calclu->DetectorType != Detector_t::Type_t::CB)
                continue;
            hCandClTime->Fill(calclu->Time,calclu->CentralElement);
            hCandClEnergy->Fill(calclu->Energy,calclu->CentralElement);
        }

}

void scratch_lheijken_checkcb::Finish()
{
}

void scratch_lheijken_checkcb::ShowResult()
{
    canvas(GetName()) << drawoption("colz")
            << hCBPhiPIDPhi
            << endc;
}

void scratch_lheijken_checkcb::CreateHistos()
{
    auto hfGeneric = new HistogramFactory("Generic",HistFac,"");
    auto hfDReadHits = new HistogramFactory("DReadHits", HistFac, "");
    auto hfClustHits = new HistogramFactory("ClustHits", HistFac, "");
    auto hfClusters = new HistogramFactory("Clusters", HistFac, "");
    auto hfCandidates = new HistogramFactory("Candidates",HistFac,"");


    const BinSettings UncalTimeBins = BinSettings(2000,-1000,1000);
    const BinSettings CalTimeBins = BinSettings(1500,-750,750);
    const BinSettings UncalEnergyBins = BinSettings(10000);
    const BinSettings CalEnergyBins = BinSettings(2000,0,4000);
    const BinSettings CalEnergyBinsS = BinSettings(1000,0,2000);
//    // logarithmic from about 10 (raw units, threshold is at 16) to 2^14, which is maximum ADC value. We need to full energy range for good TimeWalk correction!
//    const BinSettings UncalEnergyBinsLog(500,std::log10(10),std::log10(1 << 14));
    //const auto CalTimeBins2 = BinSettings::RoundToBinSize({100,-50,50}, calibration::converter::Gains::CATCH_TDC);
//    const auto CalTimeBins3 = BinSettings(100,-50,50);
    const BinSettings PhiBins(500,-180.,180.);
    double pidphiminch0 = std_ext::radian_to_degree(pid_detector->GetPosition(0).Phi()) - std_ext::radian_to_degree(pid_detector->dPhi(0))/2.;
    double pidphimin = pidphiminch0 - std_ext::radian_to_degree(pid_detector->dPhi(0))*24;
    double pidphimax = pidphiminch0+ std_ext::radian_to_degree(pid_detector->dPhi(0))*24;
    std::cout<<"phimin="<<pidphimin<<", phimax="<<pidphimax<<", dPhi="<<std_ext::radian_to_degree(pid_detector->dPhi(0))<<std::endl;
    const BinSettings PIDPhiBins(48,pidphimin,pidphimax);

    hTrigRefTiming = hfGeneric->makeTH1D("Trigger reference timing","Trigger reference timing [ns]","",CalTimeBins,"hTrigRefTiming",true);

    hDRHCBesum = hfDReadHits->makeTH1D("CB DRH esum","CB esum [MeV]","",CalEnergyBinsS,"hDRHCBesum",true);
    hDRHCBesum_1Cl = hfDReadHits->makeTH1D("CB DRH esum, only 1 cluster","CB esum [MeV]","",CalEnergyBinsS,"hDRHCBesum_1Cl",true);
    hDRHCBesum_2Cl = hfDReadHits->makeTH1D("CB DRH esum, only 2 cluster","CB esum [MeV]","",CalEnergyBinsS,"hDRHCBesum_2Cl",true);
    hDRHCBesum_3Cl = hfDReadHits->makeTH1D("CB DRH esum, only 3 cluster","CB esum [MeV]","",CalEnergyBinsS,"hDRHCBesum_3Cl",true);
    hDRHCBesum_4Cl = hfDReadHits->makeTH1D("CB DRH esum, only 4 cluster","CB esum [MeV]","",CalEnergyBinsS,"hDRHCBesum_4Cl",true);
    hDRHCBesum_5Cl = hfDReadHits->makeTH1D("CB DRH esum, only 5 cluster","CB esum [MeV]","",CalEnergyBinsS,"hDRHCBesum_5Cl",true);
    hDRHUncalTimeAll = hfDReadHits->makeTH2D("CB DRH UncalTime all", "uncalibrated time", "CB channel", UncalTimeBins, BinSettings(720),"hDRHUncalTimeAll", true);
    hDRHCalTimeAll = hfDReadHits->makeTH2D("CB DRH CalTime all", "calibrated time", "CB channel", CalTimeBins, BinSettings(720),"hDRHCalTimeAll", true);
    hDRHUncalTimeFirst = hfDReadHits->makeTH2D("CB DRH UncalTime first", "uncalibrated time", "CB channel", UncalTimeBins, BinSettings(720),"hDRHUncalTimeFirst", true);
    hDRHCalTimeFirst = hfDReadHits->makeTH2D("CB DRH CalTime first", "calibrated time", "CB channel", CalTimeBins, BinSettings(720),"hDRHCalTimeFirst", true);
    hDRHUncalEnergy = hfDReadHits->makeTH2D("CB DRH UncalEnergy", "uncalibrated energy", "CB channel", UncalEnergyBins, BinSettings(720),"hDRHUncalEnergy", true);
    hDRHCalEnergy = hfDReadHits->makeTH2D("CB DRH CalEnergy", "calibrated energy", "CB channel", CalEnergyBins, BinSettings(720),"hDRHCalEnergy", true);
    hDRHUncalTimeMult = hfDReadHits->makeTH2D("CB DRH UncalTime multiplicity","Pulse multiplicity","CB channel",BinSettings(20),BinSettings(720),"hDRHUncalTimeMult", true);
    hDRHUncalEnergyMult = hfDReadHits->makeTH2D("CB DRH UncalEnergy multiplicity","Pulse multiplicity","CB channel",BinSettings(20),BinSettings(720),"hDRHUncalEnergyMult", true);
    hDRHUnCalEn_wTiming = hfDReadHits->makeTH2D("CB DRH UncalEnergy with timing", "uncalibrated energy", "CB channel", UncalEnergyBins, BinSettings(720),"hDRHUnCalEn_wTiming", true);
    hDRHUnCalEn_woTiming = hfDReadHits->makeTH2D("CB DRH UncalEnergy without timing", "uncalibrated energy", "CB channel", UncalEnergyBins, BinSettings(720),"hDRHUnCalEn_woTiming", true);
    hDRHCalEn_wTiming = hfDReadHits->makeTH2D("CB DRH CalEnergy with timing", "calibrated energy", "CB channel", CalEnergyBins, BinSettings(720),"hDRHCalEn_wTiming", true);
    hDRHCalEn_woTiming = hfDReadHits->makeTH2D("CB DRH CalEnergy without timing", "calibrated energy", "CB channel", CalEnergyBins, BinSettings(720),"hDRHCalEn_woTiming", true);

    hCHTime = hfClustHits->makeTH2D("CB ClustHit Time","Time","CB channel",CalTimeBins,BinSettings(720),"hCHTime",true);
    hCHEnergy = hfClustHits->makeTH2D("CB ClustHit Energy","Energy","CB channel",CalEnergyBins,BinSettings(720),"hCHEnergy",true);
//    hCHTimeRawE = hfClustHits->makeTH3D("CB ClustHit Time vs RawEnergy","Time","log_{10}(RawEnergy)","Channel",CalTimeBins3,UncalEnergyBinsLog,BinSettings(720),"hCHTimeRawE",true);

    hCBPhiPIDPhi = hfClusters->makeTH2D("CB Phi vs PID Phi","CB Phi","PID Phi",PhiBins,PIDPhiBins,"hCBPhiPIDPhi",true);
    hClTime = hfClusters->makeTH2D("CB Cluster Time","Time","CB channel",CalTimeBins,BinSettings(720),"hClTime",true);
    hClEnergy = hfClusters->makeTH2D("CB Cluster Energy","Energy","CB channel",CalEnergyBins,BinSettings(720),"hClEnergy",true);
    hClTimeVsEnergy = hfClusters->makeTH2D("Cluster Time vs Energy","Cluster Time","Cluster Energy",CalTimeBins,CalEnergyBins,"hClTimeVsEnergy",true);

    hCandClTime = hfCandidates->makeTH2D("CB Candidate Time","Time","CB channel",CalTimeBins,BinSettings(720),"hCandClTime",true);
    hCandClEnergy = hfCandidates->makeTH2D("CB Candidate Energy","Energy","CB channel",CalEnergyBins,BinSettings(720),"hCandClEnergy",true);


}
AUTO_REGISTER_PHYSICS(scratch_lheijken_checkcb)

