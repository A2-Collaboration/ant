#include "FitGausPol3.h"

#include "TF1.h"
#include "TH1.h"
#include "TF1Knobs.h"

#include "base/Logger.h"


void ant::calibration::gui::FitGausPol3::sync()
{
    signal->SetParameters(&(combinded->GetParameters()[0]));
    bg->SetParameters(    &(combinded->GetParameters()[3]));
    setRange(signal,GetRange());
    setRange(bg,GetRange());
}

ant::calibration::gui::FitGausPol3::FitGausPol3()
{
    signal = new TF1("","gaus");
    signal->SetLineColor(kRed);

    bg = new TF1("","pol3");
    bg->SetLineColor(kBlue);

    combinded = new TF1("","gaus(0)+pol3(3)");
    combinded->SetLineColor(kGreen);

    SetRange(ant::interval<double>(100,250));
    combinded->SetParName(0,"A");
    combinded->SetParName(1,"x_{0}");
    combinded->SetParName(2,"#sigma");
    combinded->SetParName(3,"p_{0}");
    combinded->SetParName(4,"p_{1}");
    combinded->SetParName(5,"p_{2}");
    combinded->SetParName(6,"p_{3}");

    combinded->SetParameter(0, 1000);
    combinded->SetParameter(1, 135);
    combinded->SetParameter(2, 8);
    combinded->SetParameter(3, 1);
    combinded->SetParameter(4, 1);
    combinded->SetParameter(5, 1);
    combinded->SetParameter(6, 0.1);

    combinded->SetParLimits(1, 115, 140);
    combinded->SetParLimits(2, 5, 50);
    combinded->FixParameter(6, 0);

    sync();

    Addknob<KnobsTF1::ParameterKnob>(combinded->GetParName(0), combinded, 0, GUIElementDescription::GUI_Type::slider_horizontal);
    Addknob<KnobsTF1::ParameterKnob>(combinded->GetParName(1), combinded, 1, GUIElementDescription::GUI_Type::slider_vertical);
    Addknob<KnobsTF1::ReferenceParameterKnob>(combinded->GetParName(2), combinded, 2, 1, GUIElementDescription::GUI_Type::slider_vertical);
    Addknob<KnobsTF1::RangeKnob>("min",combinded,KnobsTF1::RangeKnob::RangeEndType::lower);
    Addknob<KnobsTF1::RangeKnob>("max",combinded,KnobsTF1::RangeKnob::RangeEndType::upper);
}

ant::calibration::gui::FitGausPol3::~FitGausPol3()
{
    delete signal;
    delete bg;
    delete combinded;
}

void ant::calibration::gui::FitGausPol3::Draw()
{
    signal->Draw("same");
    bg->Draw("same");
    combinded->Draw("same");
}

void ant::calibration::gui::FitGausPol3::Fit(TH1* hist)
{
    hist->Fit(combinded, "RBQN");
    sync();
}

void ant::calibration::gui::FitGausPol3::SetRange(ant::interval<double> i)
{
    setRange(combinded, i);
    setRange(signal, i);
    setRange(bg, i);
}

ant::interval<double> ant::calibration::gui::FitGausPol3::GetRange() const
{
    return getRange(combinded);
}

void ant::calibration::gui::FitGausPol3::Sync()
{
    sync();
}

void ant::calibration::gui::FitGausPol3::SetPoints(int n)
{
    signal->SetNpx(n);
    bg->SetNpx(n);
    combinded->SetNpx(n);
}

std::vector<double> ant::calibration::gui::FitGausPol3::Save() const
{
    SavedState_t params;
    params.reserve(2+combinded->GetNpar());
    saveTF1(combinded,params);

    return params;
}

void ant::calibration::gui::FitGausPol3::Load(const std::vector<double> &data)
{
    if(data.size() != std::size_t(2+combinded->GetNpar())) {
        LOG(WARNING) << "Can't load parametes";
        return;
    }

    SavedState_t::const_iterator p = data.begin();
    loadTF1(p, combinded);
    SetRange(getRange(combinded));
    sync();
}

