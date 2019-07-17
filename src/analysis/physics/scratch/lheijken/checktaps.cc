#include "checktaps.h"

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
 * checking info from taps
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_lheijken_checktaps::scratch_lheijken_checktaps(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get())
{
    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    taps_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
    veto_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPSVeto>();

    CreateHistos();

}


void scratch_lheijken_checktaps::ProcessEvent(const TEvent& event, manager_t&)
{

    triggersimu.ProcessEvent(event);

    const double TriggerRefTime = triggersimu.GetRefTiming();
    hTrigRefTiming->Fill(TriggerRefTime);

    //-- Loop over DetectorReadhits in event and extract the ones from TAPS to look at
    int elemhastiming[taps_detector->GetNChannels()];
    memset(elemhastiming, 0, taps_detector->GetNChannels()*sizeof(int));
    const auto& readhits = event.Reconstructed().DetectorReadHits;
    for(const TDetectorReadHit& readhit : readhits) {
        if(readhit.DetectorType != Detector_t::Type_t::TAPS)
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
            }
            hDRHUncalEnergyMult->Fill(emult,readhit.Channel);
        }
    }
    for(const TDetectorReadHit& readhit : readhits) {
        if(readhit.DetectorType != Detector_t::Type_t::TAPS)
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
    //-- Loop over Cluster and Clusterhits in event and extract the ones from TAPS to look at
    auto& clusters = event.Reconstructed().Clusters;
    for (const TCluster& cluster : clusters){
        if(cluster.DetectorType != Detector_t::Type_t::TAPS )
            continue;
        for(const TClusterHit& hit : cluster.Hits){
            hCHTime->Fill(hit.Time,hit.Channel);
            hCHEnergy->Fill(hit.Energy,hit.Channel);          
        }
        hClTime->Fill(cluster.Time,cluster.CentralElement);
        hClEnergy->Fill(cluster.Energy,cluster.CentralElement);
        hClTimeVsEnergy->Fill(cluster.Time,cluster.Energy);
    }

    //-- Loop over candidates and look at their TAPS clusters
    const auto& cands = event.Reconstructed().Candidates;
    for (const TCandidate& cand : cands){
            const auto calclu = cand.FindCaloCluster();
            if(calclu->DetectorType != Detector_t::Type_t::TAPS)
                continue;
            hCandClTime->Fill(calclu->Time,calclu->CentralElement);
            hCandClEnergy->Fill(calclu->Energy,calclu->CentralElement);
        }

}

void scratch_lheijken_checktaps::Finish()
{
}

void scratch_lheijken_checktaps::ShowResult()
{

}

void scratch_lheijken_checktaps::CreateHistos()
{
    auto hfGeneric = new HistogramFactory("Generic",HistFac,"");
    auto hfDReadHits = new HistogramFactory("DReadHits", HistFac, "");
    auto hfClustHits = new HistogramFactory("ClustHits", HistFac, "");
    auto hfClusters = new HistogramFactory("Clusters", HistFac, "");
    auto hfCandidates = new HistogramFactory("Candidates",HistFac,"");


    const BinSettings UncalTimeBins = BinSettings(2500,-1000,4000);
    const BinSettings CalTimeBins = BinSettings(1500,-750,750);
    const BinSettings UncalEnergyBins = BinSettings(10000);
    const BinSettings CalEnergyBins = BinSettings(2000,0,4000);

    hTrigRefTiming = hfGeneric->makeTH1D("Trigger reference timing","Trigger reference timing [ns]","",CalTimeBins,"hTrigRefTiming",true);

    hDRHUncalTimeAll = hfDReadHits->makeTH2D("TAPS DRH UncalTime all", "uncalibrated time", "TAPS channel", UncalTimeBins, BinSettings(720),"hDRHUncalTimeAll", true);
    hDRHCalTimeAll = hfDReadHits->makeTH2D("TAPS DRH CalTime all", "calibrated time", "TAPS channel", CalTimeBins, BinSettings(720),"hDRHCalTimeAll", true);
    hDRHUncalTimeFirst = hfDReadHits->makeTH2D("TAPS DRH UncalTime first", "uncalibrated time", "TAPS channel", UncalTimeBins, BinSettings(720),"hDRHUncalTimeFirst", true);
    hDRHCalTimeFirst = hfDReadHits->makeTH2D("TAPS DRH CalTime first", "calibrated time", "TAPS channel", CalTimeBins, BinSettings(720),"hDRHCalTimeFirst", true);
    hDRHUncalEnergy = hfDReadHits->makeTH2D("TAPS DRH UncalEnergy", "uncalibrated energy", "TAPS channel", UncalEnergyBins, BinSettings(720),"hDRHUncalEnergy", true);
    hDRHCalEnergy = hfDReadHits->makeTH2D("TAPS DRH CalEnergy", "calibrated energy", "TAPS channel", CalEnergyBins, BinSettings(720),"hDRHCalEnergy", true);
    hDRHUncalTimeMult = hfDReadHits->makeTH2D("TAPS DRH UncalTime multiplicity","Pulse multiplicity","TAPS channel",BinSettings(20),BinSettings(720),"hDRHUncalTimeMult", true);
    hDRHUncalEnergyMult = hfDReadHits->makeTH2D("TAPS DRH UncalEnergy multiplicity","Pulse multiplicity","TAPS channel",BinSettings(20),BinSettings(720),"hDRHUncalEnergyMult", true);
    hDRHUnCalEn_wTiming = hfDReadHits->makeTH2D("TAPS DRH UncalEnergy with timing", "uncalibrated energy", "TAPS channel", UncalEnergyBins, BinSettings(720),"hDRHUnCalEn_wTiming", true);
    hDRHUnCalEn_woTiming = hfDReadHits->makeTH2D("TAPS DRH UncalEnergy without timing", "uncalibrated energy", "TAPS channel", UncalEnergyBins, BinSettings(720),"hDRHUnCalEn_woTiming", true);
    hDRHCalEn_wTiming = hfDReadHits->makeTH2D("TAPS DRH CalEnergy with timing", "calibrated energy", "TAPS channel", CalEnergyBins, BinSettings(720),"hDRHCalEn_wTiming", true);
    hDRHCalEn_woTiming = hfDReadHits->makeTH2D("TAPS DRH CalEnergy without timing", "calibrated energy", "TAPS channel", CalEnergyBins, BinSettings(720),"hDRHCalEn_woTiming", true);

    hCHTime = hfClustHits->makeTH2D("TAPS ClustHit Time","Time","TAPS channel",CalTimeBins,BinSettings(720),"hCHTime",true);
    hCHEnergy = hfClustHits->makeTH2D("TAPS ClustHit Energy","Energy","TAPS channel",CalEnergyBins,BinSettings(720),"hCHEnergy",true);

    hClTime = hfClusters->makeTH2D("TAPS Cluster Time","Time","TAPS channel",CalTimeBins,BinSettings(720),"hClTime",true);
    hClEnergy = hfClusters->makeTH2D("TAPS Cluster Energy","Energy","TAPS channel",CalEnergyBins,BinSettings(720),"hClEnergy",true);
    hClTimeVsEnergy = hfClusters->makeTH2D("Cluster Time vs Energy","Cluster Time","Cluster Energy",CalTimeBins,CalEnergyBins,"hClTimeVsEnergy",true);

    hCandClTime = hfCandidates->makeTH2D("TAPS Candidate Time","Time","TAPS channel",CalTimeBins,BinSettings(720),"hCandClTime",true);
    hCandClEnergy = hfCandidates->makeTH2D("TAPS Candidate Energy","Energy","TAPS channel",CalEnergyBins,BinSettings(720),"hCandClEnergy",true);


}
AUTO_REGISTER_PHYSICS(scratch_lheijken_checktaps)

