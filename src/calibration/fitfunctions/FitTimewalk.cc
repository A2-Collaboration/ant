#include "FitTimewalk.h"

#include "base/interval.h"
#include "base/Logger.h"
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
    func->SetParName(2, "E_{0}");

    AddKnob<KnobsTF1::RangedParameterKnob>(func->GetParName(2), func, 2,
                                           KnobsTF1::RangedParameterKnob::ConstraintType::lowerThanMin,
                                           IndicatorProperties::Type_t::slider_vertical);

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
    FitFunction::doFit(hist, func);
}

void FitTimewalk::SetDefaults(TH1*)
{
    func->SetParameter(0,  -25);
    func->SetParameter(1,   55);
    func->SetParameter(2,   -5);
    func->SetParameter(3, 0.15);

    func->SetParLimits(2, -100, 0);

    SetRange({5, 400});
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

double FitTimewalk::Eval(double energy)
{
    return func->Eval(energy);
}


