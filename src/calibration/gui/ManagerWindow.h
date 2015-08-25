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

struct ManagerWindowMode {
    ManagerWindowMode() :
        gotoNextSlice(true),
        alwaysDisplayFit(false),
        channelStep(1)
    {}

    bool gotoNextSlice;
    bool alwaysDisplayFit;
    int  channelStep;
};

class ManagerWindow : public TGMainFrame
{
private:
    std::list<CalCanvas*> canvases;
    TGHorizontalFrame* frame_canvases;
    TGStatusBar* statusbar;
    std::map<EKeySym, TGTextButton*> keys;
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