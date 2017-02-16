#pragma once

#include "FitFunction.h"

namespace ant {
namespace calibration {
namespace gui {

class FitWeibullLandauPol1: public PeakingFitFunction {

    using transformation_t = std::function<double(double, TF1*)>;
    transformation_t MPV_trafo;
    transformation_t MPV_trafo_inverse;

protected:
    TF1* signal;
    TF1* bg;

public:
    FitWeibullLandauPol1();

    virtual ~FitWeibullLandauPol1();

    virtual void Draw() override;

    virtual void Fit(TH1* hist) override;
    virtual void FitSignal(TH1* hist) override;
    virtual void FitBackground(TH1* hist) override;
    virtual void SetDefaults(TH1* hist) override;

    virtual void SetRange(ant::interval<double> i) override;
    virtual ant::interval<double> GetRange() const override;

    virtual void Sync() override;

    virtual SavedState_t Save() const override;
    virtual void Load(const SavedState_t &data) override;

    virtual double GetPeakPosition() const override;
    virtual double GetPeakWidth() const override;
    virtual double SignalToBackground(const double x) const override;

    double GetWeibullPeak() const;
    double GetWeibullWidth() const;
    double GetLandauPeak() const;
    double GetLandauWidth() const;
};

}
}
}

