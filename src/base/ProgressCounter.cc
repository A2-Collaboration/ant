#include "ProgressCounter.h"

#include "std_ext/math.h"

#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace ant;
using namespace std;

long int ProgressCounter::Interval = 0;

ProgressCounter::registry_t ProgressCounter::registry;
ProgressCounter::timepoint_t ProgressCounter::last_now = ProgressCounter::clock_t::now();

void ProgressCounter::Tick()
{
    if(Interval<=0)
        return;
    const auto now = clock_t::now();
    std::chrono::duration<double> elapsed = now - last_now;
    if(elapsed.count()<Interval)
        return;
    for(auto c : registry)
        c->Updater(elapsed);
    last_now = now;
}

ProgressCounter::ProgressCounter(ProgressCounter::Updater_t updater) :
    Updater(updater)
{
    registry.push_back(this);
}

ProgressCounter::~ProgressCounter() {
    auto it = std::find(registry.begin(), registry.end(), this);
    registry.erase(it);
}

double ProgressCounter::GetTotalSecs() const {
    return static_cast<chrono::duration<double>>(clock_t::now() - first_now).count();
}

string ProgressCounter::TimeToStr(double secs)
{
    if(!isfinite(secs))
        return "??";

    const auto hours_left = std::floor(secs / 3600);
    secs -= hours_left * 3600;
    const auto mins_left = std::floor(secs / 60);
    secs -= mins_left * 60;

    std::stringstream ss_ETA;
    if(hours_left>0)
        ss_ETA << hours_left << ":";
    if(mins_left>0)
        ss_ETA << setw(2) << setfill('0') << mins_left << ":";
    ss_ETA << setw(2) << setfill('0') << std::floor(secs);

    return ss_ETA.str();
}
