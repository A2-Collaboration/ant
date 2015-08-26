#pragma once

#include "calibration/gui/Indicator_traits.h"
#include "calibration/fitfunctions/FitFunction.h"

#include <list>
#include <memory>
#include <string>
#include <stack>

#include "TCanvas.h"


class TH1;
class TRootCanvas;
class TGStatusBar;

namespace ant {
namespace calibration {
namespace gui {

struct CalCanvasMode {
    CalCanvasMode() :
        gotoNextRange(true),
        alwaysDisplayFit(false),
        channelStep(1)
    {}

    bool gotoNextRange;
    bool alwaysDisplayFit;
    int  channelStep;
};

class CalCanvas : public TCanvas, public update_notify_traits {
protected:
    Viewport getViewport();
    void ClearIndicators();
    void SetupGUI();

    Indicator* MakeVerticalIndicatorLine(IndicatorKnob& knob);
    Indicator* MakeHorizontalIndicatorLine(IndicatorKnob& knob);
    Indicator* MakeGUIElement(IndicatorKnob& knob);

    std::list<Indicator*> indicators;
    FitFunction* func = nullptr;
    TH1* hist = nullptr;

    virtual void HandleKeypress(const char key);

    std::stack<FitFunction::SavedState_t> UndoStack;

    TRootCanvas* rootcanvas = nullptr;
    CalCanvasMode* gui_mode = nullptr;
    TGStatusBar* statusbar = nullptr;


public:
    CalCanvas(const std::string& name, const std::string& title);

    CalCanvas(const std::string& name, Int_t winid) :
        TCanvas(name.c_str(), 10, 10, winid) {}

    virtual ~CalCanvas();

    virtual void Show(TH1* h, FitFunction* f);

    virtual void UpdateMe() override;


    virtual void Fit();
    void SetDefaults();
    virtual void UndoPush();
    virtual void UndoPop();

    virtual void ConnectReturnFunc(const char* receiver_class, void* receiver, const char* slot);
    virtual void LinkGUIMode(CalCanvasMode* guimode_);

    virtual void ConnectStatusBar(TGStatusBar* statusbar_);

    virtual void Execute(const char *method, const char *params, Int_t *error);

    virtual void Update() override;

    /**
     * @brief HandleInput: Override default to catch keyboard inputs
     * @param button equals kKeyPress
     * @param x the key char which was pressed
     * @param y
     */
    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;

    /**
     * @brief ShowGuidelines remove the snapping guidelines completely
     */
    // do not mark it override since older ROOT versions don't have this routine
    virtual void ShowGuidelines(TObject*, const Int_t, const char, const bool) {}

    virtual void ProcessedEvent(Int_t event, Int_t x, Int_t y, TObject *selected) override;
};
}
}
}
