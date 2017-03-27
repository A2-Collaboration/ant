#pragma once

#include "FitFunction.h"

namespace ant {
namespace calibration {
namespace gui {

class FitGausPol0: public PeakingFitFunction {

public:
    struct p {
        constexpr static auto Offset  = 3;
        constexpr static auto Height  = 0;
        constexpr static auto Pos     = 1;
        constexpr static auto Sigma   = 2;
    };

    FitGausPol0();

    virtual ~FitGausPol0();

    virtual void Draw() override;

    virtual void Fit(TH1* hist) override;
    virtual void SetDefaults(TH1* hist) override;

    virtual void SetRange(ant::interval<double> i) override;
    virtual ant::interval<double> GetRange() const override;

    virtual SavedState_t Save() const override;
    virtual void Load(const SavedState_t &data) override;

    virtual double GetPeakPosition() const override;
    virtual double GetPeakWidth() const override;

    double SignalToBackground(const double x) const override;
};

}
}
}

