#pragma once

#include "FitFunction.h"

namespace ant {
namespace calibration {
namespace gui {

class FitVetoBand: public FitFunction {

    void EnsureParameterLimits();

protected:
    TF1* signal;
    TF1* bg;

public:
    FitVetoBand();
    ~FitVetoBand();

    void Draw() override;

    void Fit(TH1* hist) override;
    void FitSignal(TH1* hist) override;
    void FitBackground(TH1* hist) override;
    void SetDefaults(TH1* hist) override;

    void SetRange(ant::interval<double> i) override;
    ant::interval<double> GetRange() const override;

    virtual void Sync() override;

    SavedState_t Save() const override;
    void Load(const SavedState_t &data) override;

    double Eval(const double energy) const;

    double EvalReference(const double energy) const;

};

}
}
}

