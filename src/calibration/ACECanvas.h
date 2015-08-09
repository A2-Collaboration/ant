#pragma once

#include "CalibrationEditor.h"

#include <string>
#include <map>
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

    void makeCalHist();

    void change_state(state_t newstate);

    void loadFile(const std::string& fileName);

    void HandleKeypress(const char key);

    void addSomeRandomData();

public:
    ACECanvas(const std::string& fileName);
//    SaveToFile(const std::string& fileName);

    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;


};

}
