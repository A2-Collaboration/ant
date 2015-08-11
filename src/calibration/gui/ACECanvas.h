#pragma once

#include "calibration/CalibrationEditor.h"
#include "base/interval.h"

#include <string>
#include <map>
#include <set>
#include <functional>
#include "Rtypes.h"
#include "analysis/plot/root_draw.h"

#include "TCanvas.h"
#include "TPaveText.h"

class TH2D;


namespace ant {

class TextPad: public root_drawable_traits
{
protected:
    double xMargin;
    double yMargin;

    TPaveText* paveText;

    // root_drawable_traits interface
public:
    TextPad(double xmargin = 0.05, double ymargin = 0.05);
    TextPad(std::list<std::string>& lines,
            double xmargin = 0.05, double ymargin = 0.05);
    virtual void Add(const std::string& line) { paveText->AddText(line.c_str());}
    void Draw(const std::string& option) const override;
};



namespace calibration {
namespace gui {




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
        reduceToValid,
        save,
        loadID,
        close
    };

    state_t             state;
    std::string         fileName;
    std::string         currentCalID;
    calibration::Editor ed;

    std::unique_ptr<ant::TextPad> display;
    TH2D*                         calHist;

    std::set<std::uint32_t> indexMemory;
    interval<std::uint32_t> indexInterVal;
    bool                    intervalStartSet;

    void update_modified();

    void makeCalHist();
    void loadFile(const std::string& fileName);

    void change_state(state_t newstate);

    void HandleKeypress(const char key);

    void addSomeRandomData();

    void updateCalHist();
    void markLine(Int_t y);
    void removeAllinIndexMemory();
    void removeInterValFromMemory();
    void markInterval(Int_t y);
    void unFillLine(uint32_t lineNumber);
    void fillLine(uint32_t lineNumber);
    void expandAllinIndexMemory();
    void markUnValid();

    void saveToFile(const std::string& fileName);

    void loadCalibration();
    void changeCalibrationID();
public:
    ACECanvas(const std::string& fileName);
    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;
};

}
}
}
