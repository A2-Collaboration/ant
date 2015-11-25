#include "ProgressCounter.h"

#include "base/Logger.h"
#include <sstream>
#include <string>
#include <iomanip>

using namespace ant;
using namespace std;

using duration_t = chrono::duration<double>;

ProgressCounter::ProgressCounter(double start, double stop):
    last_output(chrono::system_clock::now()), x(start), max(stop), last_x(x)
{}

void ProgressCounter::SetInterval(double i)
{
    t = i;
}

string TimeToStr(unsigned secs) {

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

bool ProgressCounter::Update(double pos)
{
    x = pos;

    const auto now = chrono::system_clock::now();
    const duration_t elapsed_seconds = now - last_output;

    if(elapsed_seconds.count() >= t) {
        percent_sec  = (x - last_x) / elapsed_seconds.count();
        sec_emaining = (max - x) / percent_sec;
        last_output = now;
        last_x = x;
        return true;
    }
    return false;
}

double ProgressCounter::PercentDone() const
{
    return double(x)/double(max)*100.0;
}

double ProgressCounter::SecondsLeft() const
{
    return sec_emaining;
}


ostream&operator<<(ostream& stream, const ProgressCounter& counter)
{
    stream << setw(2) << std::setprecision(4) << counter.PercentDone() << " % done, ETA: " << TimeToStr(unsigned(counter.SecondsLeft()));
    return stream;
}
