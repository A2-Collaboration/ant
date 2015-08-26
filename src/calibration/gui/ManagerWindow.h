#pragma once

#include "Manager_traits.h"

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

class ManagerWindow : public TGMainFrame, public ManagerWindow_traits
{
private:
    std::list<CalCanvas*> canvases;
    TGHorizontalFrame* frame_canvases;
    TGStatusBar* statusbar;
    std::map<EKeySym, TGTextButton*> keys;
    Manager* manager;
    void CreateToolbar(TGVerticalFrame* frame);
    void UpdateLayout();
public:
    ManagerWindow(Manager* manager_);
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~ManagerWindow();
    virtual gui::CalCanvas* AddCalCanvas(const std::string& name = "") override;

    ManagerWindowMode Mode;

    ManagerWindow(const ManagerWindow&) = delete;
    ManagerWindow& operator=(const ManagerWindow&) = delete;
};
}
}
}