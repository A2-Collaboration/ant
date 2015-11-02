#include "FitGausPol3.h"


#include "BaseFunctions.h"

#include "base/Logger.h"

#include "base/TF1Ext.h"

#include "TF1.h"
#include "TH1.h"

using namespace std;
using namespace ant::calibration;

void ant::calibration::gui::FitGausPol3::sync()
{
    signal->SetParameters(&(combined->GetParameters()[0]));
    bg->SetParameters(    &(combined->GetParameters()[3]));
    setRange(signal,GetRange());
    setRange(bg,GetRange());
}

ant::calibration::gui::FitGausPol3::FitGausPol3()
{
    signal = functions::gaus::getTF1();
    signal->SetLineColor(kRed);

    bg = functions::pol<3>::getTF1();
    bg->SetLineColor(kBlue);

    combined = functions::GausPol<3>::getTF1();
    combined->SetLineColor(kGreen);

    SetRange(ant::interval<double>(100,250));
    combined->SetParName(0,"A");
    combined->SetParLimits(0, 0.0, 1E+12);
    combined->SetParName(1,"x_{0}");
    combined->SetParName(2,"#sigma");
    combined->SetParLimits(2, 0.0, 1E+12);
    combined->SetParName(3,"p_{0}");
    combined->SetParName(4,"p_{1}");
    combined->SetParName(5,"p_{2}");
    combined->SetParName(6,"p_{3}");

    signal->SetNpx(1000);
    bg->SetNpx(1000);
    combined->SetNpx(1000);

    AddKnob<KnobsTF1::ParameterKnob>(combined->GetParName(0), combined, 0, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::ParameterKnob>(combined->GetParName(1), combined, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::ReferenceParameterKnob>(combined->GetParName(2), combined, 2, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::RangeKnob>("min",combined,KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("max",combined,KnobsTF1::RangeKnob::RangeEndType::upper);
}

ant::calibration::gui::FitGausPol3::~FitGausPol3()
{
    delete signal;
    delete bg;
    delete combined;
}

void ant::calibration::gui::FitGausPol3::Draw()
{
    signal->Draw("same");
    bg->Draw("same");
    combined->Draw("same");
}

void ant::calibration::gui::FitGausPol3::Fit(TH1* hist)
{
    FitFunction::doFit(hist, combined);
    sync();
}

void gui::FitGausPol3::FitBackground(TH1* hist)
{
    const auto fixedPars = {0,1,2};
    FixParameters(combined, fixedPars);
    FitFunction::doFit(hist, combined);
    sync();
    UnFixParameters(combined, fixedPars);
}

void gui::FitGausPol3::FitSignal(TH1* hist)
{
    const auto fixedPars = {3,4,5,6};
    FixParameters(combined, fixedPars);
    FitFunction::doFit(hist, combined);
    sync();
    UnFixParameters(combined, fixedPars);
}

void ant::calibration::gui::FitGausPol3::SetDefaults(TH1 *hist)
{
    // defaults for taps baf2
    combined->SetParameter(0, hist->GetMaximum());
    combined->SetParameter(1, 135);
    combined->SetParameter(2, 8);
    combined->SetParameter(3, 1);
    combined->SetParameter(4, 1);
    combined->SetParameter(5, 1);
    combined->SetParameter(6, 0.1);

    combined->SetParLimits(1, 115, 140);
    combined->SetParLimits(2, 5, 50);

    sync();
}

void ant::calibration::gui::FitGausPol3::SetRange(ant::interval<double> i)
{
    setRange(combined, i);
    setRange(signal, i);
    setRange(bg, i);
    // x_0 peak position must be in range
    combined->SetParLimits(1, i.Start(), i.Stop());
}

ant::interval<double> ant::calibration::gui::FitGausPol3::GetRange() const
{
    return getRange(combined);
}

void ant::calibration::gui::FitGausPol3::Sync()
{
    sync();
}

std::vector<double> ant::calibration::gui::FitGausPol3::Save() const
{
    SavedState_t params;
    params.reserve(2+combined->GetNpar());
    saveTF1(combined,params);

    return params;
}

void ant::calibration::gui::FitGausPol3::Load(const std::vector<double> &data)
{
    if(data.size() != std::size_t(2+combined->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }

    SavedState_t::const_iterator p = data.begin();
    loadTF1(p, combined);
    SetRange(getRange(combined));
    sync();
}

double ant::calibration::gui::FitGausPol3::GetPeakPosition() const
{
    return combined->GetParameter(1);
}

double gui::FitGausPol3::GetPeakWidth() const
{
    return combined->GetParameter(2);
}

