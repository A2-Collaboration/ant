#pragma once

#include "calibration/gui/manager/GUIInterface.h"
#include <memory>

class TH1;
class TObject;

namespace ant {
namespace calibration {
namespace gui {

class FitGausPol3;
class CalCanvas;

class DebugModule: public GUIClientInterface {
protected:
    std::shared_ptr<FitGausPol3> func;

public:
    DebugModule ();
    virtual ~DebugModule();

    std::string GetHistogramName() const;
    FitStatus Fit(CalCanvas* c, TH1* hist, unsigned channel);
    void StoreResult(unsigned channel);
    FitStatus Finish();
    unsigned GetNumberOfChannels();
};
}
}
}
