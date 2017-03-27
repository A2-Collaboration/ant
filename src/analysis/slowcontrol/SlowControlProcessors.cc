#include "SlowControlProcessors.h"

using namespace std;
using namespace ant::analysis::slowcontrol;

#define DEFINE_PROCESSOR(proc) \
    const shared_ptr<processor::proc> Processors::proc = make_shared<processor::proc>();

DEFINE_PROCESSOR(EPT_Scalers)
DEFINE_PROCESSOR(EPT_Or)
DEFINE_PROCESSOR(Tagger_Scalers)
DEFINE_PROCESSOR(Tagger_Or)
DEFINE_PROCESSOR(Beampolmon)
DEFINE_PROCESSOR(ExpTrigger)
DEFINE_PROCESSOR(Beam)


