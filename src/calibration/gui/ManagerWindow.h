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

class ManagerWindow : public TGMainFrame
{
private:
    std::list<CalCanvas*> canvases;
    TGHorizontalFrame* frame_canvases;
    TGStatusBar* statusbar;
    std::map<EKeySym, TGTextButton*> keys;
public:
    ManagerWindow(const TGWindow*p, UInt_t w, UInt_t h);
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~ManagerWindow();
    gui::CalCanvas* AddCalCanvas(const std::string& name = "");
};
}
}
}