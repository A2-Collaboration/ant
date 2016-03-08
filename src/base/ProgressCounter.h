#pragma once

#include "std_ext/memory.h"


#include <chrono>
#include <thread>
#include <functional>
#include <string>

namespace ant {

struct ProgressCounter {

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
    void work();

    Updater_t Updater;
    std::unique_ptr<std::thread> worker;

    using clock_t = std::chrono::steady_clock;
    std::chrono::time_point<clock_t> first_now = clock_t::now();
    std::chrono::time_point<clock_t> last_now = clock_t::now();

};



}