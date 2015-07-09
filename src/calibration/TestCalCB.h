#ifndef TESTCALCB_H
#define TESTCALCB_H

#include "BaseCalModule.h"

class TH1;

namespace ant {
namespace calibration {

class TestCalCB: public BaseCalibrationModule {
protected:
    TH1* ggIM = nullptr;

public:
    TestCalCB();
    virtual ~TestCalCB() = default;

    // CalibrationApply_traits interface
public:
    void ApplyTo();

    // Physics interface
public:
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}

#endif
