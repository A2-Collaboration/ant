#include "Slowcontrol.h"
#include "base/std_ext/memory.h"

#include "tree/TSlowControl.h"

#include <map>
#include <stdexcept>
#include <limits>

using namespace std;
using namespace ant;
using namespace ant::analysis::data;

bool SlowcontrolRequestable::isRequested() const
{
    return requested;
}

void SlowcontrolRequestable::Request()
{
    requested = true;
}



template <> ant::analysis::data::ReqestableVariable<double>::ReqestableVariable():
    data(std::numeric_limits<double>::quiet_NaN())
{

}

template <> ant::analysis::data::ReqestableVariable<std::int64_t>::ReqestableVariable():
    data(std::numeric_limits<int64_t>::lowest())
{

}
