#include "SlowControlVariables.h"

using namespace std;
using namespace ant::analysis::slowcontrol;

// ensure the list is created before the variables
list<VariablePtr> Variables::All;

#define DEFINE_VARIABLE(var) \
    const shared_ptr<variable::var> Variables::var = make_shared<variable::var>(); \
    Variables::AddToAll addToAll ## var (Variables::var);

DEFINE_VARIABLE(TaggerScalers)
DEFINE_VARIABLE(PhotonFlux)
