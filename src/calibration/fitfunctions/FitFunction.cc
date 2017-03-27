#include "FitFunction.h"
#include "BaseFunctions.h"


#include "base/interval.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;


ant::interval<double> FitFunction::getRange(const TF1* func)
{
    interval<double> i(std_ext::NaN, std_ext::NaN);
    func->GetRange(i.Start(), i.Stop());
    return i;
}

void FitFunction::setRange(TF1* func, const ant::interval<double>& i)
{
    func->SetRange(i.Start(), i.Stop());
}

void FitFunction::doFit(TH1* hist, TF1* func)
{
    hist->Fit(func,"RBNQW");
}

void FitFunction::saveTF1(const TF1 *func, SavedState_t &out)
{
    auto range = getRange(func);
    out.push_back(range.Start());
    out.push_back(range.Stop());

    for(int i=0; i <func->GetNpar(); ++i ) {
        out.push_back(func->GetParameter(i));
    }
}

void FitFunction::loadTF1(SavedState_t::const_iterator &data_pos, TF1 *func)
{
    setRange(func,{*(data_pos++),*(data_pos++)});

    std::copy(data_pos,data_pos+func->GetNpar(),func->GetParameters());
}

FitFunction::~FitFunction()
{}

double FitFunction::Chi2NDF() const
{
    return Chi2() / NDF();
}

double FitFunction::Chi2() const
{
    return func->GetChisquare();
}

double FitFunction::NDF() const
{
    return func->GetNDF();
}

PeakingFitFunction::PeakingFitFunction(): FitFunction(){}

double PeakingFitFunction::SignalToBackground(const double x) const
{
    return std::abs(func->Eval(x));
}

bool PeakingFitFunction::EndsMatch(const double relative_epsilon) const
{
    const auto range = GetRange();
    return      std::abs(SignalToBackground(range.Start())) < relative_epsilon
            &&  std::abs(SignalToBackground(range.Stop()) ) < relative_epsilon;
}
