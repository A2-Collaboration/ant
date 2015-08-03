#pragma once

#include "calibration/gui/manager/GUIInterface.h"
#include <memory>

class TH1;
class TObject;


namespace ant {
namespace calibration {
namespace gui {

class FitGausPol3;

class DebugModule: public GUIClientInrerface {
protected:
    TQObject* guiinstance = nullptr;
    std::shared_ptr<FitGausPol3> func;

public:
    DebugModule ();
    virtual ~DebugModule();

    std::string GetHistogramName() const;
    FitStatus Fit(TH1* hist, unsigned channel);
    void StoreResult(unsigned channel);
    FitStatus Finish();
    TQObject* GetGUIInstance();
    unsigned GetNumberOfChannels();
};
}
}
}
