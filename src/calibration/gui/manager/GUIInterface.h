#pragma once
#include <string>

class TH1;
class TQObject;

namespace ant {
namespace calibration {
namespace gui {

//class GUIMasterInterface {
//public:
//    virtual
//};

class GUIClientInterface {
public:
    enum class FitStatus {
        FitOK,
        GUIWait
    };

    virtual std::string GetHistogramName() const =0;
    virtual FitStatus Fit(TH1* hist, unsigned channel) =0;
    virtual void StoreResult(unsigned channel) =0;
    virtual FitStatus Finish() =0;
    virtual TQObject* GetGUIInstance() =0;
    virtual unsigned GetNumberOfChannels() =0;
};

}
}
}
