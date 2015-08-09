#pragma once

#include "CalibrationEditor.h"
#include "base/interval.h"

#include <string>
#include <map>
#include <set>
#include <functional>
#include "Rtypes.h"

#include "TCanvas.h"

class TH2D;


namespace ant {


class ACECanvas: public TCanvas
{
private:
    enum class state_t
    {
        init,
        base,
        remove,
        cut,
        expand,
        close
    };

    state_t state;
    std::string fileName;
    std::string currentCalID;
    calibration::Editor ed;
    TH2D* calHist;

    std::set<std::uint32_t> stepMemory;

    interval<std::uint32_t> stepInterVal;
    bool intervalStartSet;

    void makeCalHist();

    void change_state(state_t newstate);

    void loadFile(const std::string& fileName);

    void HandleKeypress(const char key);

    void addSomeRandomData();

    void updateCalHist();
    void markLine(Int_t y);
    void removeAllinStepMemory();
    void removeInterValFromMemory();
    void markInterval(Int_t y);
public:
    ACECanvas(const std::string& fileName);
//    SaveToFile(const std::string& fileName);

    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;


};

}
