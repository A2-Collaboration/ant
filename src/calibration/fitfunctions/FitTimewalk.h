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

class FitTimewalk: public FitFunction {
protected:
    TF1* func = nullptr;

public:
    FitTimewalk();

    virtual ~FitTimewalk();

    virtual void Draw() override;

    virtual void Fit(TH1* hist) override;
    virtual void SetDefaults(TH1* hist) override;

    virtual void SetRange(ant::interval<double> i) override;
    virtual ant::interval<double> GetRange() const override;

    virtual SavedState_t Save() const override;
    virtual void Load(const SavedState_t &data) override;
};

}
}
}

