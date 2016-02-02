#include "SlowControl.h"
#include "base/std_ext/memory.h"

#include "tree/TSlowControl.h"

#include <map>
#include <stdexcept>
#include <limits>

using namespace std;
using namespace ant::analysis::slowcontrol;

bool SlowControlRequestable::isRequested() const
{
    return requested;
}

void SlowControlRequestable::Request()
{
    requested = true;
}



template <> ReqestableVariable<double>::ReqestableVariable():
    data(std::numeric_limits<double>::quiet_NaN())
{

}

template <> ReqestableVariable<std::int64_t>::ReqestableVariable():
    data(std::numeric_limits<int64_t>::lowest())
{

}
