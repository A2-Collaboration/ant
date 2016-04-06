#pragma once

#include <base/interval.h>
#include <cmath>

namespace ant {

template <typename T>
T step(const ant::interval<T> i, const unsigned n, const unsigned p) {
    return i.Start() + i.Length() / (n-1) * std::min(n-1,p);
}

struct stepper {
    const ant::interval<double>& i;
    const unsigned n;
    unsigned p=0;
    double value;

    stepper(const interval<double>& interv, const unsigned steps): i(interv), n(steps), value(step(i,n,0)) {}

    bool Next() {

        if(p<n) {
            ++p;
            value = step(i,n,p);
            return true;
        } else
            return false;
    }

    bool Done() const {
        return p >= n;
    }

};

}
