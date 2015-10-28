#pragma once

#include "calibration/gui/Indicator_traits.h"
#include "calibration/fitfunctions/FitFunction.h"

#include <list>
#include <memory>
#include <string>
#include <stack>
#include <limits>

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

    TGStatusBar* statusbar = nullptr;

    struct axis_settings_t {
        double x_min = std::numeric_limits<double>::quiet_NaN();
        double x_max = std::numeric_limits<double>::quiet_NaN();
        double y_min = std::numeric_limits<double>::quiet_NaN();
        double y_max = std::numeric_limits<double>::quiet_NaN();
        bool HaveX() const {
            return std::isfinite(x_min) && std::isfinite(x_max);
        }
        bool HaveY() const {
            return std::isfinite(y_min) && std::isfinite(y_max);
        }
    };
    axis_settings_t axis_settings;

public:
    // constructor for simple canvas
    CalCanvas(const std::string& name) :
        TCanvas(name.c_str()) {}


    // constructor for embedded canvas
    CalCanvas(const std::string& name, Int_t winid) :
        TCanvas(name.c_str(), 10, 10, winid) {}

    virtual ~CalCanvas();

    virtual void Show(TH1* h, FitFunction* f, bool preserveYaxis = false);

    virtual void UpdateMe() override;

    enum class FitType_t {
        Total,
        Signal,
        Background
    };

    virtual void Fit(const FitType_t type=FitType_t::Total);
    void SetDefaults();
    virtual void UndoPush();
    virtual void UndoPop();


    virtual void Update() override;
    void Update(bool fromhist);

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
