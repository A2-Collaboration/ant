#include "Slowcontrol.h"
#include "base/std_ext/memory.h"

#include "tree/TSlowControl.h"

#include <map>
#include <stdexcept>

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

