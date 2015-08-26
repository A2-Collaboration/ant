#pragma once

#include "TGFrame.h"
#include "KeySymbols.h"

#include <functional>
#include <list>
#include <map>

class TGStatusBar;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvas;
class Manager;

struct ManagerWindowMode {
    ManagerWindowMode() :
        gotoNextSlice(true),
        autoContinue(true),
        showEachFit(true),
        channelStep(1)
    {}

    bool gotoNextSlice;
    bool autoContinue;
    bool showEachFit;
    int  channelStep;
};

class ManagerWindow : public TGMainFrame
{
private:
    std::list<CalCanvas*> canvases;
    TGHorizontalFrame* frame_canvases;
    TGStatusBar* statusbar;
    std::map<EKeySym, TGTextButton*> keys;
    Manager* manager;
    void CreateToolbar(TGVerticalFrame* frame);
public:
    ManagerWindow(const TGWindow* p, UInt_t w, UInt_t h);
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~ManagerWindow();
    gui::CalCanvas* AddCalCanvas(const std::string& name = "");

    ManagerWindowMode Mode;

    ManagerWindow(const ManagerWindow&) = delete;
    ManagerWindow& operator=(const ManagerWindow&) = delete;
};
}
}
}