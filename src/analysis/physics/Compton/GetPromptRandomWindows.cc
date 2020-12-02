#include "GetPromptRandomWindows.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

// Defining the contructor for the Compton class
GetPromptRandomWindows::GetPromptRandomWindows(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    // Histograms are created here but not filled

    const BinSettings time_bins(1000, -100, 100);

    h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
                                    "t [ns]","#",     // xlabel, ylabel
                                    time_bins,        // binnings
                                    "h_TaggerTime"    // ROOT object name, auto-generated if omitted
                                    );
    h_TaggerTimeCBSubtraction = HistFac.makeTH1D("Tagger Time with CB Time Subtracted",
                                    "t [ns]","#",
                                    time_bins,
                                    "h_TaggerTimeCBSubtraction"
                                    );
}

void GetPromptRandomWindows::ProcessEvent(const TEvent& event, manager_t&)
{
    // Runs the ProcessEvent function for TriggerSimulation
    // does the calculations
    triggersimu.ProcessEvent(event);

    for (const auto& taggerhit : event.Reconstructed().TaggerHits)
    {
        // Apply trigger simulation to tagger hits
        // This subtracts a weighted time from the CB (see wiki)
        const auto& CorrectedTaggerTime =
                triggersimu.GetCorrectedTaggerTime(taggerhit);

        promptrandom.SetTaggerTime(CorrectedTaggerTime);
        h_TaggerTime->Fill(taggerhit.Time);
        h_TaggerTimeCBSubtraction->Fill(CorrectedTaggerTime);
    }
}

void GetPromptRandomWindows::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time Plot")
            << h_TaggerTime                      // Tagger hit times
            << h_TaggerTimeCBSubtraction         // Tagger hit times with the subtraction
                                                 // -> used to identify PR windows
            << endc; // draws the canvas
}

// A macro that registers the GetPromptRandomWindows class with Ant so that
// you can call this class using "Ant -p GetPromptRandomWindows"
AUTO_REGISTER_PHYSICS(GetPromptRandomWindows)
