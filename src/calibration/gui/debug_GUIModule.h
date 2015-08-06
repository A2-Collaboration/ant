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

    std::string GetHistogramName() const override;
    FitStatus Fit(CalCanvas* c, TH1* hist, unsigned channel) override;
    void StoreResult(unsigned channel) override;
    FitStatus Finish(CalCanvas* c) override;
    void StoreFinish() override;
    unsigned GetNumberOfChannels();
};
}
}
}
