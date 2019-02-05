#include "checktagger.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "slowcontrol/SlowControlVariables.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "base/WrapTFile.h"

#include <iostream>
#include <memory>

/**
 * checking info from tagger
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_lheijken_checktagger::scratch_lheijken_checktagger(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get())
{
    CreateHistos();
    TaggHitTree.CreateBranches(HistFac.makeTTree("tagghit_variables"));

    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();

    slowcontrol::Variables::Beam->Request();
}


void scratch_lheijken_checktagger::ProcessEvent(const TEvent& event, manager_t&)
{

    triggersimu.ProcessEvent(event);

    const double TriggerRefTime = triggersimu.GetRefTiming();
    hTriggerRefTiming->Fill(TriggerRefTime);

    //-- Loop over ReadHits and fetch the tagger multiplicities from there
    int taggmult[368];
    int firstchannel = 500;
    for(int i=0; i<368; i++) taggmult[i]=-1;
    for(const TDetectorReadHit& dethit : event.Reconstructed().DetectorReadHits) {
        if(dethit.DetectorType != Detector_t::Type_t::Tagger)
            continue;
        if(dethit.ChannelType != Channel_t::Type_t::Timing)
            continue;
        if(dethit.Values.size()<1)
            continue;
        hTimeMultiplicity->Fill(dethit.Values.size(), dethit.Channel);
        taggmult[dethit.Channel] = dethit.Values.size();
        if(firstchannel > (int)dethit.Channel)
            firstchannel = dethit.Channel;
    }

    //-- Loop over tagger hits
    //--- Some counters to keep track of tagger multiplicity
    int itiming = 0;
    int previouschannel = firstchannel;
    //--- Some counters and vectors to store all times and count #hits in the first TDC module
    vector <double> HitTimes;
    vector <double> HitChannel;
    int nrHitsInTDCMod0=0;
    int nrHitsInTDCMod1=0;
    int nrHitsInTDCMod2=0;
    for (const auto& tHit: event.Reconstructed().TaggerHits) {
        if(previouschannel == (int)tHit.Channel)
            itiming++;
        else{
            hTimeMultiplicityFromTHits->Fill(itiming, previouschannel);
            hTimeMultTaggVsReadHits->Fill(itiming,taggmult[previouschannel]);
            itiming=1;
        }

        //---- Fill histograms
        hTime->Fill(tHit.Time, tHit.Channel);
        hTimeToTriggerRef->Fill(tHit.Time - TriggerRefTime, tHit.Channel);
        //---- Fill tree
        TaggHitTree.EventNumber = event.Reconstructed().ID.Lower;
        TaggHitTree.Time = tHit.Time;
        TaggHitTree.Channel = tHit.Channel;
        TaggHitTree.ChannelMult = taggmult[tHit.Channel];
        TaggHitTree.TimingInd = itiming;
        TaggHitTree.Tree->Fill();

        previouschannel = tHit.Channel;

        //---- Check which tagger channel belongs to which TDC module
        hTestChanNrMod[(unsigned)tagger_detector->GetTDCSector(tHit.Channel)]->Fill(tHit.Channel);

        //---- Store the time and the channel number
        HitTimes.push_back(tHit.Time);
        HitChannel.push_back(tHit.Channel);
        if((unsigned)tagger_detector->GetTDCSector(tHit.Channel) < 1)
            nrHitsInTDCMod0++;
        else if((unsigned)tagger_detector->GetTDCSector(tHit.Channel) < 2)
            nrHitsInTDCMod1++;
        else
            nrHitsInTDCMod2++;
    }
    hTimeMultiplicityFromTHits->Fill(itiming, previouschannel);
    hTimeMultTaggVsReadHits->Fill(itiming,taggmult[previouschannel]);
    itiming=0;

    //--- Fill the check on nr hits in each TDC module
    hNrHitInTDCModPerEv[0]->Fill(nrHitsInTDCMod0);
    hNrHitInTDCModPerEv[1]->Fill(nrHitsInTDCMod1);
    hNrHitInTDCModPerEv[2]->Fill(nrHitsInTDCMod2);
    for(int i=0; i<((int)HitTimes.size());i++){
        if(nrHitsInTDCMod0<170)
            hTime_cutNrHitInTDCMod0[0]->Fill(HitTimes[i],HitChannel[i]);
        if(nrHitsInTDCMod0<160)
            hTime_cutNrHitInTDCMod0[1]->Fill(HitTimes[i],HitChannel[i]);
        if(nrHitsInTDCMod0<140)
            hTime_cutNrHitInTDCMod0[2]->Fill(HitTimes[i],HitChannel[i]);

      }

    if(slowcontrol::Variables::Beam->HasChanged())
    {
        hFaradayCupScaler->Fill(slowcontrol::Variables::Beam->GetFaradayCup());
    }
}

void scratch_lheijken_checktagger::Finish()
{
}

void scratch_lheijken_checktagger::ShowResult()
{

}

void scratch_lheijken_checktagger::CreateHistos()
{
    const BinSettings TimeBins = BinSettings(1500,-750,750);

    hTime = HistFac.makeTH2D("Tagger Time", "time [ns]", "Tagger channel", TimeBins, BinSettings(368),"hTime");
    hTimeToTriggerRef = HistFac.makeTH2D("Tagger Time relative to TriggerRef", "time", "Tagger channel",TimeBins, BinSettings(368),"hTimeToTriggerRef");
    hTriggerRefTiming = HistFac.makeTH1D("CB - Trigger timing", "time [ns]", "#", BinSettings(500,-50,50), "hTriggerRefTiming");
    hTimeMultiplicity = HistFac.makeTH2D("Tagger Hit Multiplicity", "multiplicity", "Tagger channel", BinSettings(20), BinSettings(368), "hTimeMultiplicity");
    hTimeMultiplicityFromTHits = HistFac.makeTH2D("Tagger Hit Multiplicity from hits", "multiplicity", "Tagger channel", BinSettings(20), BinSettings(368), "hTimeMultiplicityFromTHits");
    hTimeMultTaggVsReadHits = HistFac.makeTH2D("Tagger Hit Multiplicity, tagghits vs readhits","mult tagghits","mult readhits",BinSettings(20),BinSettings(20),"hTimeMultTaggVsReadHits");
    hNrHitInTDCModPerEv[0] = HistFac.makeTH1D("Nr hits in TDCmodule 0 per event","Nr hits in TDCmodule 0 per event","#",BinSettings(300,-0.5,299.5),"hNrHitInTDCMod0PerEv");
    hNrHitInTDCModPerEv[1] = HistFac.makeTH1D("Nr hits in TDCmodule 1 per event","Nr hits in TDCmodule 1 per event","#",BinSettings(300,-0.5,299.5),"hNrHitInTDCMod1PerEv");
    hNrHitInTDCModPerEv[2] = HistFac.makeTH1D("Nr hits in TDCmodule 2 per event","Nr hits in TDCmodule 2 per event","#",BinSettings(300,-0.5,299.5),"hNrHitInTDCMod2PerEv");
    hTime_cutNrHitInTDCMod0[0] = HistFac.makeTH2D("Tagger Time w cut on nr hits in TDCmod 0/ev < 170", "time [ns]", "Tagger channel", TimeBins, BinSettings(368),"hTime_cutNrHitInTDCMod0_1");
    hTime_cutNrHitInTDCMod0[1] = HistFac.makeTH2D("Tagger Time w cut on nr hits in TDCmod 0/ev < 160", "time [ns]", "Tagger channel", TimeBins, BinSettings(368),"hTime_cutNrHitInTDCMod0_2");
    hTime_cutNrHitInTDCMod0[2] = HistFac.makeTH2D("Tagger Time w cut on nr hits in TDCmod 0/ev < 140", "time [ns]", "Tagger channel", TimeBins, BinSettings(368),"hTime_cutNrHitInTDCMod0_3");

    hTestChanNrMod[0] = HistFac.makeTH1D("Channel nr module 1","Tagger channel","#",BinSettings(368),"hTestChanNrMod_1");
    hTestChanNrMod[1] = HistFac.makeTH1D("Channel nr module 2","Tagger channel","#",BinSettings(368),"hTestChanNrMod_2");
    hTestChanNrMod[2] = HistFac.makeTH1D("Channel nr module 3","Tagger channel","#",BinSettings(368),"hTestChanNrMod_3");

    hFaradayCupScaler = HistFac.makeTH1D("FaradayCupScaler","Beam current [pA]","#",BinSettings(3000,0.,30000),"hFaradayCupScaler");

}
AUTO_REGISTER_PHYSICS(scratch_lheijken_checktagger)

