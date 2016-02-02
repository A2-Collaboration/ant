#include "EventFilter.h"

#include "TH1D.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

void SetBinLabel(TH1* hist, const int bin, const std::string& label) {
    hist->GetXaxis()->SetBinLabel(bin, label.c_str());
}

EventFilter::EventFilter(const string& name, OptionsPtr opts):
    Physics(name, opts),
    nCands(opts->Get<interval<size_t>>("nCands", {0, 100})),
    CBEsum(opts->Get<double>("CBEsum",    550.0))
{
    steps = HistFac.makeTH1D("Filter Steps", "", "", BinSettings(3,0,3));
    SetBinLabel(steps,1, formatter() << "Input");
    SetBinLabel(steps,2, formatter() << "CBEsum > " << CBEsum);
    SetBinLabel(steps,3, formatter() << "nCands in " << nCands);

    VLOG(3) << "CBEsum > " << CBEsum;
    VLOG(3) << "nCands in " << nCands;
}

EventFilter::~EventFilter()
{}

void EventFilter::ProcessEvent(const TEvent& event, Physics::manager_t& manager)
{
    const auto& data = *(event.Reconstructed);

    steps->SetBinContent(1, steps->GetBinContent(1)+1);

    if(data.Trigger.CBEnergySum < CBEsum)
        return;

    steps->SetBinContent(2, steps->GetBinContent(2)+1);

    if(!nCands.Contains(data.Candidates.size()))
        return;

    steps->SetBinContent(3, steps->GetBinContent(3)+1);

    manager.SaveEvent();
}

AUTO_REGISTER_PHYSICS(EventFilter)
