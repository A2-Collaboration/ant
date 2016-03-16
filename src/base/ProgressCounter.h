#pragma once

#include "std_ext/memory.h"


#include <chrono>
#include <functional>
#include <string>
#include <list>

namespace ant {

struct ProgressCounter {

    static void Tick();

    using Updater_t = std::function<void(std::chrono::duration<double>)>;
    ProgressCounter(Updater_t updater);
    ~ProgressCounter();

    double GetTotalSecs() const;

    /**
     * @brief TimeToStr converts seconds to HH:MM:SS string
     * @param secs seconds to convert
     * @return string with format HH:MM:SS
     */
    static std::string TimeToStr(double secs);

    /**
     * @brief Interval for progress output in seconds
     */
    static long int Interval;

protected:
    Updater_t Updater;

    using clock_t = std::chrono::steady_clock;
    using timepoint_t = std::chrono::time_point<clock_t>;
    timepoint_t first_now = clock_t::now();

    using registry_t = std::list<const ProgressCounter*>;
    static registry_t  registry;
    static timepoint_t last_now;
};



}