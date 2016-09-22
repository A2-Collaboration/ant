#include "base/vec/vec2.h"

namespace ant {
namespace math {

struct LineFct {
    double m;
    double b;

    LineFct(const vec2& p1, const vec2& p2) :
        m((p2.y - p1.y) / (p2.x - p1.x)),
        b(p1.y - m * p1.x)
    {}

    double operator() (const double& x) const noexcept { return m*x+b; }
};

}}
