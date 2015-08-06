#pragma once

#include "Manager_traits.h"
#include <memory>

class TH1;
class TObject;

namespace ant {
namespace calibration {
namespace gui {

class FitGausPol3;
class CalCanvas;

class DebugModule : public gui::Manager_traits {
protected:
    std::shared_ptr<FitGausPol3> func;
    CalCanvas* canvas  = nullptr;
    TH1* projection = nullptr;

public:
    DebugModule ();
    virtual ~DebugModule();

    std::string GetHistogramName() const override;
    unsigned GetNumberOfChannels() const override;
    std::list<CalCanvas*> GetCanvases() const override;
    void InitGUI() override;

    void StartRange(const interval<TID>& range) override;

    bool DoFit(TH1* hist, unsigned channel) override;
    void DisplayFit() override;
    void StoreFit(unsigned channel) override;

    bool FinishRange() override;
    void StoreFinishRange(const interval<TID>& range) override;
};

}}} // ant::calibration::gui
