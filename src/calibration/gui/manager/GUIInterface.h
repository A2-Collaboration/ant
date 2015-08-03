#pragma once
#include <string>

class TH1;

namespace ant {
namespace calibration {
namespace gui {

//class GUIMasterInterface {
//public:
//    virtual
//};

class GUIClientInrerface {
public:
    enum class FitStatus {
        FitOK,
        GUIWait
    };

    virtual std::string GetHistogramName() const =0;
    virtual FitStatus Fit(TH1* hist) =0;
    virtual FitStatus Finish() =0;
};

}
}
}
