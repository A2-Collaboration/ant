#pragma once

#include "FitFunction.h"

namespace ant {
namespace calibration {
namespace gui {

class FitGausPol3: public PeakingFitFunction {
protected:
    TF1* signal;
    TF1* bg;
    TF1* combined;

    void sync();


public:
    FitGausPol3();
    virtual ~FitGausPol3();

    void Draw();
    void Fit(TH1* hist);
    virtual void SetDefaults(TH1 *hist) override;
    void SetRange(ant::interval<double> i);
    ant::interval<double> GetRange() const;
    virtual void Sync() override;

    std::vector<double> Save() const;
    void Load(const std::vector<double> &data);

    virtual double GetPeakPosition() const override;

};

}
}
}
