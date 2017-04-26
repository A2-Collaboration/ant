#pragma once

#include "base/std_ext/math.h"
#include "expconfig/ExpConfig.h"

#include <random>

namespace ant {

struct TEvent;
struct TTaggerHit;

namespace analysis {
namespace utils {

class TriggerSimulation {

    struct info_t {
        bool   hasTriggered;
        double CBEnergySum;
        double CBTiming;
        info_t() {
            Reset(); // ensure reset on startup
        }
        void Reset() {
            hasTriggered = false;
            CBEnergySum = std_ext::NaN;
            CBTiming = std_ext::NaN;
        }
        bool IsSane() {
            return std::isfinite(CBEnergySum) &&
                    std::isfinite(CBTiming);
        }
    };

    info_t info;

    using config_t = expconfig::Setup_traits::triggersimu_config_t;
    const config_t config;

    std::unique_ptr<std::default_random_engine> random_gen;
    std::normal_distribution<double> random_CBESum_threshold;

public:

    TriggerSimulation();

    /**
     * @brief ProcessEvent inspects the full event to tell trigger decision
     * @param event the event under investigation
     * @return true if successful, false on failure
     */
    bool ProcessEvent(const TEvent& event);

    /**
     * @brief HasTriggered returns true if the experiment would have accepted this event
     * @return the trigger decision
     */
    bool   HasTriggered() const { return info.hasTriggered; }

    /**
     * @brief GetCBEsum tries to calculate the energy sum
     * @return Crystal Ball energy sum in MeV
     */
    double GetCBEnergySum() const { return info.CBEnergySum; }

    /**
     * @brief GetRefTiming tries to calculate the reference timing,
     * typically given by the analog energy sum
     * @return Crystal Ball timing offset in nanoseconds
     */
    double GetRefTiming() const { return info.CBTiming; }

    /**
     * @brief GetCorrectedTaggerTime improves the timing resolution by removing the trigger offset
     * @param taggerhit the taggerhit to be corrected
     * @return corrected tagger timing in ns
     * @note Typically calculated as "Taggertime - RefTiming", but might differ for various beamtimes
     */
    double GetCorrectedTaggerTime(const TTaggerHit& taggerhit) const;
};

}
}
}