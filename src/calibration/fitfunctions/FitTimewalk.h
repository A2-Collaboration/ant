#pragma once

#include "FitFunction.h"

namespace ant {
namespace calibration {
namespace gui {

class FitTimewalk: public FitFunction {

    void EnsureParameterLimits();
    bool loaded = false;

public:
    FitTimewalk();
    ~FitTimewalk();

    void Draw() override;

    void Fit(TH1* hist) override;
    void FitSignal(TH1* hist) override;
    void FitBackground(TH1* hist) override;
    void SetDefaults(TH1* hist) override;

    void SetRange(ant::interval<double> i) override;
    ant::interval<double> GetRange() const override;

    SavedState_t Save() const override;
    void Load(const SavedState_t &data) override;

    double Eval(double energy);

};

}
}
}

