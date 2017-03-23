#pragma once

#include "ManagerWindow_traits.h"

#include "TGFrame.h"
#include "KeySymbols.h"

#include <functional>
#include <list>
#include <map>

class TGStatusBar;
class TGWidget;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvas;
class Manager;
class ProgressBar;

class ManagerWindow : public TGMainFrame,
        public ManagerWindowGUI_traits
{
private:
    std::list<CalCanvas*> canvases;
    TGHorizontalFrame* frame_canvases = nullptr;
    TGHorizontalFrame* frame_extraflags = nullptr;

    TGStatusBar* statusbar = nullptr;
    ProgressBar* progress_channel = nullptr;
    ProgressBar* progress_slice = nullptr;

    std::list<TGWidget*> nonfinish_widgets;

    std::map<EKeySym, TGTextButton*> keys;
    Manager& manager;


    void CreateToolbar(TGVerticalFrame* frame);
    void UpdateLayout();
    bool running = false;
    void RunManager();

    Mode_t mode;

public:
    ManagerWindow(Manager& manager_);
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~ManagerWindow();

    gui::CalCanvas* AddCalCanvas(const std::string& name = "") override;
    void AddCheckBox(const std::string &label, bool& flag) override;
    void AddNumberEntry(const std::string &label, double& number) override;
    void AddNumberEntry(const std::string& label, double initial_number, std::function<void(const TGNumberEntry&)> callback) override;

    Mode_t& GetMode() override { return mode; }
    void SetProgressMax(unsigned slices, unsigned channels) override;
    void SetProgress(unsigned slice, unsigned channel) override;
    void SetFinishMode(bool flag) override;

    ManagerWindow(const ManagerWindow&) = delete;
    ManagerWindow& operator=(const ManagerWindow&) = delete;
};
}
}
}
