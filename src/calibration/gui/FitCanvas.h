#pragma once

#include "calibration/gui/GUIbase.h"

#include <list>
#include <memory>
#include <string>

#include "TCanvas.h"


class TH1;


namespace ant {
namespace calibration {
namespace gui {

class FitFunction;

class CalCanvas : public TCanvas, public update_notify_traits {
protected:
    Viewport getViewport();
    void ClearInidators();
    void SetupGUI();

    GUIIndicator* MakeVerticalIndicatorLine(VirtualKnob& knob);
    GUIIndicator* MakeHorizontalIndicatorLine(VirtualKnob& knob);
    GUIIndicator* MakeGUIElement(VirtualKnob& knob);

    std::list<GUIIndicator*> indicators;
    std::shared_ptr<ant::calibration::gui::FitFunction> func = nullptr;
    TH1* hist = nullptr;

public:
    CalCanvas(const std::string& name);
    virtual ~CalCanvas();

    virtual void Show(TH1* h, std::shared_ptr<FitFunction> f);

    virtual void update_me() override;

    virtual void ShowGuidelines(TObject*, const Int_t, const char, const bool) override;

    virtual void Fit();

    virtual void Execute(const char *method, const char *params, Int_t *error);

    virtual void Update() override;

};
}
}
}
