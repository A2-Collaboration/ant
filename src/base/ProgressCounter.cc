#include "ProgressCounter.h"

#include "base/Logger.h"
#include "std_ext/memory.h"
#include "std_ext/math.h"

#include <sstream>
#include <iomanip>
#include <chrono>
#include <mutex>
#include <map>

using namespace ant;
using namespace std;

long int ProgressCounter::Interval = 0;

// remember the running threads
// in this map, protect it by mutex
std::mutex running_mutex;
std::map<const thread*, bool> running;

ProgressCounter::ProgressCounter(ProgressCounter::Updater_t updater) :
    Updater(updater)
{
    // silently ignore progress output
    if(Interval<=0)
        return;

    std::lock_guard<std::mutex> lock(running_mutex);
    worker = std_ext::make_unique<std::thread>([this] () {
        while(true) {
            std::this_thread::sleep_for(std::chrono::seconds(Interval));
            std::lock_guard<std::mutex> lock(running_mutex);
            if(!worker || !running[worker.get()])
                break;
            this->work();
        }
    });
    running[worker.get()] = true;
}

ProgressCounter::~ProgressCounter() {
    if(!worker)
        return;
    worker->detach();
    std::lock_guard<std::mutex> lock(running_mutex);
    running[worker.get()] = false;
}

double ProgressCounter::GetTotalSecs() const {
    return static_cast<chrono::duration<double>>(clock_t::now() - first_now).count();
}

string ProgressCounter::TimeToStr(double secs)
{
    const unsigned hours_left = secs / 3600;
    secs -= hours_left * 3600;
    const unsigned mins_left = secs / 60;
    secs -= mins_left * 60;

    std::stringstream ss_ETA;
    if(hours_left>0)
        ss_ETA << hours_left << ":";
    if(mins_left>0)
        ss_ETA << setw(2) << setfill('0') << mins_left << ":";
    ss_ETA << setw(2) << setfill('0') << static_cast<unsigned>(secs);

    return ss_ETA.str();
}

void ProgressCounter::work() {
    const auto now = clock_t::now();
    Updater(now - last_now);
    last_now = now;
}
