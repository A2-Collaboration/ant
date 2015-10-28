#pragma once

#include "KnobsTF1.h"
#include "calibration/gui/Indicator_traits.h"

#include "base/interval.h"
#include "base/std_ext/memory.h"
#include "base/interval.h"

#include "Rtypes.h"

#include <memory>
#include <list>
#include <string>


class TH1;
class TF1;

namespace ant {
namespace calibration {
namespace gui {


class FitFunction {
public:
    using knoblist_t = std::list<std::unique_ptr<IndicatorKnob>>;
    using SavedState_t = std::vector<double>;

protected:
    knoblist_t knobs;

    template <typename T, typename ... Args_t>
    void AddKnob(Args_t&& ... args) {
        knobs.emplace_back(std_ext::make_unique<T>(std::forward<Args_t>(args)...));
    }

    static ant::interval<double> getRange(const TF1* func);
    static void setRange(TF1* func, const ant::interval<double>& i);
    static void doFit(TH1* hist, TF1* func, size_t repeat = 2);

    static void saveTF1(const TF1* func, SavedState_t& out);
    static void loadTF1(SavedState_t::const_iterator& data_pos, TF1* func);

public:
    virtual ~FitFunction();
    virtual void Draw() =0;
    knoblist_t& GetKnobs() { return knobs; }
    virtual void Fit(TH1* hist) =0;
    virtual void FitSignal(TH1*) {}
    virtual void FitBackground(TH1*) {}

    /**
     * @brief Set/Calcualte default parameter values. The hist that will be fitted later is given to allow adaptions
     * @param hist The hist to fit later
     */
    virtual void SetDefaults(TH1* hist) =0;

    virtual void SetRange(ant::interval<double> i) =0;
    virtual ant::interval<double> GetRange() const =0;
    virtual void Sync() {}

    virtual SavedState_t Save() const =0;
    virtual void Load(const std::vector<double>& data) =0;

};

class PeakingFitFunction: public FitFunction
{
public:
    PeakingFitFunction();
    virtual double GetPeakPosition() const =0;
    virtual double GetPeakWidth() const =0;
};


}
}
}

