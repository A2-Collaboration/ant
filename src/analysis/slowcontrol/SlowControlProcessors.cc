#include "SlowControlProcessors.h"

using namespace std;
using namespace ant::analysis::slowcontrol;

#define DEFINE_PROCESSOR(proc) \
    const shared_ptr<processor::proc> Processors::proc = make_shared<processor::proc>();

DEFINE_PROCESSOR(EPT_Scalers)
DEFINE_PROCESSOR(Beampolmon)