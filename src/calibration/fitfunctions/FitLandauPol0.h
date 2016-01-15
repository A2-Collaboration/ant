#pragma once

#include "FitFunction.h"

namespace ant {
namespace calibration {
namespace gui {

class FitLandauPol0: public PeakingFitFunction {

    using transformation_t = std::function<double(double, TF1*)>;
    transformation_t MPV_trafo;
    transformation_t MPV_trafo_inverse;

public:
    FitLandauPol0();

    virtual ~FitLandauPol0();

    virtual void Draw() override;

    virtual void Fit(TH1* hist) override;
    virtual void SetDefaults(TH1* hist) override;

    virtual void SetRange(ant::interval<double> i) override;
    virtual ant::interval<double> GetRange() const override;

    virtual SavedState_t Save() const override;
    virtual void Load(const SavedState_t &data) override;

    virtual double GetPeakPosition() const override;
    virtual double GetPeakWidth() const override;
};

}
}
}

