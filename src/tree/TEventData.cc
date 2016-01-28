#include "TEventData.h"
#include "TClass.h"

using namespace std;
using namespace ant;

ostream& TEventData::Print(ostream& s) const {
    return s << "TEventData ID=" << ID;
}

void TEventData::Streamer(TBuffer& R__b) {
    if (R__b.IsReading()) {
        TEventData::Class()->ReadBuffer(R__b, this);
    } else {
        TEventData::Class()->WriteBuffer(R__b, this);
    }
}

