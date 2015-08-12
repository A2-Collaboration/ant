#pragma once

#include "FitFunction.h"

#include "base/interval.h"
#include "calibration/gui/GUIbase.h"

#include <memory>

#include "TF1Knobs.h"

#include <list>
#include <string>
#include "base/std_ext.h"
#include "base/interval.h"

#include "Rtypes.h"

class TH1;
class TF1;

namespace ant {
namespace calibration {
namespace gui {

class FitFunctionGaus: public FitFunction {
protected:
    TF1* func = nullptr;

    class MyWKnob: public VirtualKnob {
    protected:
        TF1* func = nullptr;
    public:

        MyWKnob(const std::string& n, TF1* Func);

        virtual double get() const override;
        virtual void set(double a) override;

    };

public:
    FitFunctionGaus();

    virtual ~FitFunctionGaus();

    virtual void Draw() override;

    virtual void Fit(TH1* hist) override;
    virtual void SetDefaults(TH1* hist) override;

    virtual void SetRange(ant::interval<double> i) override;
    virtual ant::interval<double> GetRange() const override;
    virtual void SetPoints(int n) override;

    virtual SavedState_t Save() const override;
    virtual void Load(const SavedState_t &data) override;

    virtual double GetPeakPosition() const override;
};

}
}
}

