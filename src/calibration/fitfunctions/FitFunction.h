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

    std::string AdditionalFitArgs;

protected:
    TF1* func = nullptr;
    knoblist_t knobs;

    template <typename T, typename ... Args_t>
    void AddKnob(Args_t&& ... args) {
        knobs.emplace_back(std_ext::make_unique<T>(std::forward<Args_t>(args)...));
    }

    static ant::interval<double> getRange(const TF1* func);
    static void setRange(TF1* func, const ant::interval<double>& i);
    void doFit(TH1* hist);

    static void saveTF1(const TF1* func, SavedState_t& out);
    static void loadTF1(SavedState_t::const_iterator& data_pos, TF1* func);

public:
    virtual ~FitFunction();
    virtual void Draw() =0;
    knoblist_t& GetKnobs() { return knobs; }
    virtual void Fit(TH1* hist) =0;
    virtual void FitSignal(TH1*) {}
    virtual void FitBackground(TH1*) {}
    void SetAdditionalFitArgs(const std::string& args) { AdditionalFitArgs = args; }

    /**
     * @brief Set/Calculate default parameter values. The hist that will be fitted later is given to allow adaptions
     * @param hist The hist to fit later
     */
    virtual void SetDefaults(TH1* hist) =0;

    virtual void SetRange(ant::interval<double> i) =0;
    virtual ant::interval<double> GetRange() const =0;
    virtual void Sync() {}

    /**
     * @brief Save the current fit parameters to a vector. Can then later be loaded again using Load()
     * @return vector containing all parameters. Internal format (meaning of each double and the size of the vector) is up to the implementation of each function.
     */
    virtual SavedState_t Save() const =0;

    /**
     * @brief Load fit parameters from a vector.
     *   Useful to load previously used ones
     * @param data vector containing the values. Internal format (meaning of each double and the size of the vector) is up to the implementation of each function.
     * @see Save()
     */
    virtual void Load(const std::vector<double>& data) =0;

    /**
     * @brief Get the reduced chi^2 (=chi^2/ndf) of last fit
     * @return
     */
    virtual double Chi2NDF() const;

    /**
     * @brief Get the Chi^2 of last fit
     * @return chi^2
     */
    virtual double Chi2() const;

    /**
     * @brief Get the Number of degrees of freedom of last fit
     * @return ndf
     */
    virtual double NDF() const;

};

/**
 * @brief Fit functions describing peak shapes, like Gaus or Lorentz with optional background function
 */
class PeakingFitFunction: public FitFunction
{
public:
    PeakingFitFunction();

    /**
     * @brief Get the Position of the peak
     * @return x_max
     */
    virtual double GetPeakPosition() const =0;

    /**
     * @brief Get the Peak Width
     * @return width
     */
    virtual double GetPeakWidth() const =0;

    /**
     * @brief Signal To Background
     * @param x Position to evaluate
     * @return (Signal + Background) / (Signal - Background)
     */
    virtual double SignalToBackground(const double x) const;

    /**
     * @brief Check if Background and Total function have the same value at the range borders
     * @param relative_epsilon maximum allowed relative difference (total/background)
     * @return true if functions are equal withing limits at the range borders
     */
    virtual bool   EndsMatch(const double relative_epsilon) const;
};


}
}
}

