#include "TEventData.h"

using namespace std;
using namespace ant;

ostream& TEventData::Print(ostream& s) const {
    return s << "TEventData ID=" << ID;
}

