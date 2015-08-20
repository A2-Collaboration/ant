#include "FitTimewalk.h"

#include "base/interval.h"
#include "base/Logger.h"
#include "TF1Knobs.h"
#include "BaseFunctions.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

FitTimewalk::FitTimewalk()
{
    func = functions::timewalk::getTF1();
    func->SetNpx(1000);
    AddKnob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitTimewalk::~FitTimewalk()
{
}

void FitTimewalk::Draw()
{
    func->Draw("same");
}

void FitTimewalk::Fit(TH1 *hist)
{
    hist->Fit(func,"RBQN");
}

void FitTimewalk::SetDefaults(TH1 *hist)
{
    if(hist) {
        func->SetParameter(0,hist->GetMaximum());
    } else {
        func->SetParameter(0,100);
    }
    double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    func->SetParameter(1,max_pos);
    func->SetParameter(2,20);
    SetRange({max_pos-60, max_pos+60});
}

void FitTimewalk::SetRange(ant::interval<double> i)
{
    setRange(func, i);
}

ant::interval<double> FitTimewalk::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitTimewalk::Save() const
{
    std::vector<double> params;

    saveTF1(func,params);

    return params;
}

void FitTimewalk::Load(const SavedState_t &data)
{
    if(data.size() != std::size_t(2+func->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }
    SavedState_t::const_iterator pos = data.begin();
    loadTF1(pos, func);

    sync();
}

double FitTimewalk::GetPeakPosition() const
{
    return func->GetParameter(1);
}

