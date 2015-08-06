#pragma once
#include <string>

class TH1;
class TQObject;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvas;

class GUIClientInterface {
public:
    enum class FitStatus {
        FitOK,
        GUIWait
    };

    virtual std::string GetHistogramName() const =0;
    virtual FitStatus Fit(CalCanvas* c, TH1* hist, unsigned channel) =0;
    virtual void StoreResult(unsigned channel) =0;
    virtual FitStatus Finish(CalCanvas* c) =0;
    virtual void StoreFinish() =0;
    virtual unsigned GetNumberOfChannels() =0;
};

}
}
}
