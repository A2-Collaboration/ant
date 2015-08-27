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
class TAxis;

namespace ant {
namespace calibration {
namespace gui {

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

    std::stack<FitFunction::SavedState_t> UndoStack;

    TRootCanvas* rootcanvas = nullptr;
    TGStatusBar* statusbar = nullptr;

    void PreserveAxis(TAxis* axis1, TAxis* axis2);

public:
    CalCanvas(const std::string& name, const std::string& title);

    CalCanvas(const std::string& name, Int_t winid);

    virtual ~CalCanvas();

    virtual void Show(TH1* h, FitFunction* f, bool preserveYaxis = false);

    virtual void UpdateMe() override;


    virtual void Fit();
    void SetDefaults();
    virtual void UndoPush();
    virtual void UndoPop();


    virtual void Update() override;



    /**
     * @brief ShowGuidelines remove the snapping guidelines completely
     */
    // do not mark it override since older ROOT versions don't have this routine
    virtual void ShowGuidelines(TObject*, const Int_t, const char, const bool) {}

    virtual void ConnectStatusBar(TGStatusBar* statusbar_);
    virtual void ProcessedEvent(Int_t event, Int_t x, Int_t y, TObject *selected) override;
};
}
}
}
