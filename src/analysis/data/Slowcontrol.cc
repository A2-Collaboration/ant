#include "Slowcontrol.h"
using namespace std;
using namespace ant::analysis::data;

bool SlowcontrolRequestable::isRequested() const
{
    return requested;
}

void SlowcontrolRequestable::Request()
{
    requested = true;
}
