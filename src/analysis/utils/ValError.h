#pragma once

#include "base/std_ext/math.h"
#include <ostream>

namespace ant {
namespace analysis {
namespace utils {

struct ValError {
    double v;
    double e;

    ValError(const double& value=std_ext::NaN, const double& error=std_ext::NaN):v(value),e(error) {}

    ValError(const ValError&) = default;
    ValError(ValError&&) = default;
    ValError& operator=(const ValError&) = default;
    ValError& operator=(ValError&&) = default;

    friend std::ostream& operator<<(std::ostream& s, const ValError& v);

    ValError operator/ (const double d) const {
        return {v/d, e/d};
    }

    ValError operator/ (const ValError& o) const {
        return {
            v / o.v,
            sqrt(std_ext::sqr(e / o.v) + std_ext::sqr(v/std_ext::sqr(o.v)*o.e))
        };
    }

    static ValError Statistical(const double v) {
        return {v, sqrt(v)};
    }
};

std::ostream& operator<<(std::ostream& s, const ValError& v) {
    s << v.v << " ± " << v.e;
    return s;
}


}}} // namespace ant::analysis::utils
