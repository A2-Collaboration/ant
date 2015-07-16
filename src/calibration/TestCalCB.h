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
    virtual void ApplyTo(const std::map< Detector_t::Type_t, std::list< TDetectorReadHit* > >& hits) override;
    virtual void ApplyTo(std::unique_ptr<TEvent>&) override {}


    // Physics interface
public:
    void ProcessEvent(const Event &event) override;
    void Finish() override;
    void ShowResult() override;

    // CalibrationUpdate_traits interface
private:
    void BuildRanges(std::list<TID> &ranges) override;
    void Update(const TID &id) override;
};

}
}

#endif
