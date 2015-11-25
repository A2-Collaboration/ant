#pragma once

#include <chrono>
#include <ostream>

namespace ant {

class ProgressCounter {
protected:
    std::chrono::time_point<std::chrono::system_clock> last_output;
    double t = 5.0; //sec

    double x = 0.0;
    const double max = 1.0;

    double sec_emaining = 0.0;
    double percent_sec = 0.0;
    double last_x = 0.0;

public:
    ProgressCounter(double start=0.0, double stop=1.0);
    void SetInterval(double i);

    bool Update(double pos);

    double PercentDone() const;
    double SecondsLeft() const;

};

}

std::ostream& operator<<(std::ostream& stream, const ant::ProgressCounter& counter);
