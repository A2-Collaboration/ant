#pragma once

#include "FitFunction.h"

namespace ant {
namespace calibration {
namespace gui {

class FitTimewalk: public FitFunction {

    void EnsureParameterLimits();

public:
    FitTimewalk();

    virtual ~FitTimewalk();

    virtual void Draw() override;

    void Fit(TH1* hist) override;
    void FitSignal(TH1* hist) override;
    void FitBackground(TH1* hist) override;
    virtual void SetDefaults(TH1* hist) override;

    virtual void SetRange(ant::interval<double> i) override;
    virtual ant::interval<double> GetRange() const override;

    virtual SavedState_t Save() const override;
    virtual void Load(const SavedState_t &data) override;

    double Eval(double energy);

};

}
}
}

