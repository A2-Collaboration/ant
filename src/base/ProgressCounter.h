#pragma once

#include "base/printable.h"

#include <chrono>
#include <ostream>



namespace ant {

class ProgressCounter : public printable_traits {
protected:
    std::chrono::time_point<std::chrono::system_clock> last_output;
    double t = 5.0; //sec

    double x = 0.0;
    const double max = 1.0;

    double sec_remaining = 0.0;
    double percent_sec = 0.0;
    double last_x = 0.0;

public:
    ProgressCounter(double start=0.0, double stop=1.0);
    void SetInterval(double i);

    bool Update(double pos);

    double PercentDone() const;
    double SecondsLeft() const;

    virtual std::ostream& Print(std::ostream& stream) const override;
};

}

