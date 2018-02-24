#include "checktagger.h"

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
}


void scratch_lheijken_checktagger::ProcessEvent(const TEvent& event, manager_t&)
{

    triggersimu.ProcessEvent(event);

    const double TriggerRefTime = triggersimu.GetRefTiming();
    hTriggerRefTiming->Fill(TriggerRefTime);

    //-- Loop over ReadHits and fetch the tagger multiplicities from there
    int taggmult[352];
    int firstchannel = 500;
    for(int i=0; i<352; i++) taggmult[i]=-1;
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
        TaggHitTree.Time = tHit.Time;
        TaggHitTree.Channel = tHit.Channel;
        TaggHitTree.ChannelMult = taggmult[tHit.Channel];
        TaggHitTree.TimingInd = itiming;
        TaggHitTree.Tree->Fill();

        previouschannel = tHit.Channel;
    }
    hTimeMultiplicityFromTHits->Fill(itiming, previouschannel);
    hTimeMultTaggVsReadHits->Fill(itiming,taggmult[previouschannel]);
    itiming=0;

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

    hTime = HistFac.makeTH2D("Tagger Time", "time [ns]", "Tagger channel", TimeBins, BinSettings(352),"hTime");
    hTimeToTriggerRef = HistFac.makeTH2D("Tagger Time relative to TriggerRef", "time", "Tagger channel",TimeBins, BinSettings(352),"hTimeToTriggerRef");
    hTriggerRefTiming = HistFac.makeTH1D("CB - Trigger timing", "time [ns]", "#", BinSettings(500,-50,50), "hTriggerRefTiming");
    hTimeMultiplicity = HistFac.makeTH2D("Tagger Hit Multiplicity", "multiplicity", "Tagger channel", BinSettings(20), BinSettings(352), "hTimeMultiplicity");
    hTimeMultiplicityFromTHits = HistFac.makeTH2D("Tagger Hit Multiplicity from hits", "multiplicity", "Tagger channel", BinSettings(20), BinSettings(352), "hTimeMultiplicityFromTHits");
    hTimeMultTaggVsReadHits = HistFac.makeTH2D("Tagger Hit Multiplicity, tagghits vs readhits","mult tagghits","mult readhits",BinSettings(20),BinSettings(20),"hTimeMultTaggVsReadHits");

}
AUTO_REGISTER_PHYSICS(scratch_lheijken_checktagger)

