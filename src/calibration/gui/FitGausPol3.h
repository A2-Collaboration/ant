#pragma once

#include "FitFunction.h"

class TF1;

namespace ant {
namespace calibration {
namespace gui {

class FitGausPol3: public FitFunction {
protected:
    TF1* signal;
    TF1* bg;
    TF1* combinded;

    void sync();


public:
    FitGausPol3();
    virtual ~FitGausPol3();

    void Draw();
    void Fit(TH1* hist);
    void SetRange(ant::interval<double> i);
    ant::interval<double> GetRange() const;
    virtual void Sync() override;
    virtual void SetPoints(int n) override;

    std::vector<double> Save() const;
    void Load(const std::vector<double> &data);
};

}
}
}
