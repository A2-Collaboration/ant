#include "SlowControlVariables.h"

using namespace std;
using namespace ant::analysis::slowcontrol;

// ensure the list is created before the variables
list<VariablePtr> Variables::All;

namespace ant {
namespace analysis {
namespace slowcontrol {
struct AddToAll {
    AddToAll(const std::shared_ptr<const Variable>& var) {
        Variables::All.push_back(std::const_pointer_cast<Variable>(var));
    }
};
}}}


#define DEFINE_VARIABLE(var) \
    const shared_ptr<const variable::var> Variables::var = make_shared<variable::var>(); \
    AddToAll addToAll ## var (Variables::var);

DEFINE_VARIABLE(TaggerScalers)
DEFINE_VARIABLE(Clocks);
DEFINE_VARIABLE(Beam);
DEFINE_VARIABLE(Pairspec);
DEFINE_VARIABLE(Trigger);
DEFINE_VARIABLE(TaggEff);
