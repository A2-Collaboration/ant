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
mutex running_mutex;
vector<bool> running;

ProgressCounter::ProgressCounter(ProgressCounter::Updater_t updater) :
    Updater(updater)
{
    // silently ignore progress output
    if(Interval<=0)
        return;

    lock_guard<mutex> lock(running_mutex);

    auto it_flag = find(running.begin(), running.end(), false);
    if(it_flag == running.end())
        it_flag = running.insert(it_flag, true);
    else
        *it_flag = true;

    // n needs to be locally copied to lambda
    auto n = distance(running.begin(), it_flag);
    worker = std_ext::make_unique<std::thread>([this, n] () {
        while(true) {
            this_thread::sleep_for(std::chrono::seconds(Interval));
            lock_guard<mutex> lock(running_mutex);
            if(!running[n])
                break;
            this->work();
        }
    });
    // remember the index
    worker_n = n;
}

ProgressCounter::~ProgressCounter() {
    if(!worker)
        return;
    worker->detach();
    std::lock_guard<std::mutex> lock(running_mutex);
    running[worker_n] = false;
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

void ProgressCounter::work() {
    const auto now = clock_t::now();
    Updater(now - last_now);
    last_now = now;
}
