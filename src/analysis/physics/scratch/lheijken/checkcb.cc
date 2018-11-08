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
    CreateHistos();

    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();

}


void scratch_lheijken_checkcb::ProcessEvent(const TEvent& event, manager_t&)
{

    triggersimu.ProcessEvent(event);

    const double TriggerRefTime = triggersimu.GetRefTiming();

    //-- Loop over DetectorReadhits in event and extract the ones from CB to look at
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

    //-- Loop over Cluster and Clusterhits in event and extract the ones from CB to look at
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
    }

}

void scratch_lheijken_checkcb::Finish()
{
}

void scratch_lheijken_checkcb::ShowResult()
{
    canvas(GetName()) << drawoption("colz")
            << hDRHUncalTimeAll
            << hDRHCalTimeAll
            << hDRHUncalEnergy
            << hDRHCalEnergy
            << endc;
}

void scratch_lheijken_checkcb::CreateHistos()
{
    auto hfDReadHits = new HistogramFactory("DReadHits", HistFac, "");
    auto hfClustHits = new HistogramFactory("ClustHits", HistFac, "");


    const BinSettings UncalTimeBins = BinSettings(2000,-1000,1000);
    const BinSettings CalTimeBins = BinSettings(1500,-750,750);
    const BinSettings UncalEnergyBins = BinSettings(10000);
    const BinSettings CalEnergyBins = BinSettings(2000,0,4000);
//    // logarithmic from about 10 (raw units, threshold is at 16) to 2^14, which is maximum ADC value. We need to full energy range for good TimeWalk correction!
//    const BinSettings UncalEnergyBinsLog(500,std::log10(10),std::log10(1 << 14));
    //const auto CalTimeBins2 = BinSettings::RoundToBinSize({100,-50,50}, calibration::converter::Gains::CATCH_TDC);
//    const auto CalTimeBins3 = BinSettings(100,-50,50);

    hDRHUncalTimeAll = hfDReadHits->makeTH2D("CB DRH UncalTime all", "uncalibrated time", "CB channel", UncalTimeBins, BinSettings(720),"hDRHUncalTimeAll", true);
    hDRHCalTimeAll = hfDReadHits->makeTH2D("CB DRH CalTime all", "calibrated time", "CB channel", CalTimeBins, BinSettings(720),"hDRHCalTimeAll", true);
    hDRHUncalTimeFirst = hfDReadHits->makeTH2D("CB DRH UncalTime first", "uncalibrated time", "CB channel", UncalTimeBins, BinSettings(720),"hDRHUncalTimeFirst", true);
    hDRHCalTimeFirst = hfDReadHits->makeTH2D("CB DRH CalTime first", "calibrated time", "CB channel", CalTimeBins, BinSettings(720),"hDRHCalTimeFirst", true);
    hDRHUncalEnergy = hfDReadHits->makeTH2D("CB DRH UncalEnergy", "uncalibrated energy", "CB channel", UncalEnergyBins, BinSettings(720),"hDRHUncalEnergy", true);
    hDRHCalEnergy = hfDReadHits->makeTH2D("CB DRH CalEnergy", "calibrated energy", "CB channel", CalEnergyBins, BinSettings(720),"hDRHCalEnergy", true);
    hDRHUncalTimeMult = hfDReadHits->makeTH2D("CB DRH UncalTime multiplicity","Pulse multiplicity","CB channel",BinSettings(20),BinSettings(720),"hDRHUncalTimeMult", true);
    hDRHUncalEnergyMult = hfDReadHits->makeTH2D("CB DRH UncalEnergy multiplicity","Pulse multiplicity","CB channel",BinSettings(20),BinSettings(720),"hDRHUncalEnergyMult", true);

    hCHTime = hfClustHits->makeTH2D("CB ClustHit Time","Time","CB channel",CalTimeBins,BinSettings(720),"hCHTime",true);
    hCHEnergy = hfClustHits->makeTH2D("CB ClustHit Energy","Energy","CB channel",CalEnergyBins,BinSettings(720),"hCHEnergy",true);
//    hCHTimeRawE = hfClustHits->makeTH3D("CB ClustHit Time vs RawEnergy","Time","log_{10}(RawEnergy)","Channel",CalTimeBins3,UncalEnergyBinsLog,BinSettings(720),"hCHTimeRawE",true);

}
AUTO_REGISTER_PHYSICS(scratch_lheijken_checkcb)

