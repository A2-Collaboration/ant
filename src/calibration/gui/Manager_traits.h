#pragma once

#include "base/interval.h"
#include "tree/THeaderInfo.h"

#include <string>
#include <list>

class TH1;
class TQObject;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvas;


class ManagerWindow_traits {
public:
    virtual gui::CalCanvas* AddCalCanvas(const std::string& name = "") =0;
    virtual void HideCalCanvas(gui::CalCanvas* canvas) =0;
};

class Manager_traits {
private:
    const std::string name;
public:
    Manager_traits(const std::string& name_) :
        name(name_) {}
    virtual ~Manager_traits() {}
    virtual std::string GetName() const { return name; }

    virtual std::string GetHistogramName() const =0;
    virtual unsigned GetNumberOfChannels() const =0;

    virtual void InitCanvases(gui::ManagerWindow_traits* window) =0;

    virtual void StartRange(const interval<TID>& range) =0;

    enum class DoFitReturn_t {
        Next, Display, Skip
    };
    virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) =0;
    virtual void DisplayFit() =0;
    virtual void StoreFit(unsigned channel) =0;

    virtual bool FinishRange() =0;
    virtual void StoreFinishRange(const interval<TID>& range) =0;
};


}
}
}
