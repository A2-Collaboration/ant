#include "Cluster.h"
#include <cmath>
#include "base/std_ext/math.h"

using namespace ant::analysis::data;
using namespace ant::std_ext;

double Cluster::GetPSARadius() const
{
    return sqrt(sqr(Energy) + sqr(ShortEnergy));
}

double Cluster::GetPSAAngle() const
{
    return atan2(ShortEnergy, Energy);
}
