#include "TEvent.h"
#include "TClass.h"

#include "base/cereal/types/polymorphic.hpp"
#include "base/cereal/types/memory.hpp"
#include "base/cereal/types/vector.hpp"
#include "base/cereal/types/list.hpp"                 // hidden in TParticleTree_t ...
#include "base/cereal/archives/portable_binary.hpp"


#include "base/Logger.h"

#include <sstream>

#include <iostream>

using namespace std;
using namespace ant;

void TEvent::Streamer(TBuffer& R__b) {

    stringstream ss;

    if (R__b.IsReading()) {

        string s;
        R__b.ReadStdString(s);
        VLOG(9) << "Read  " << s.length() << " bytes" << endl;
        ss << s;
        cereal::PortableBinaryInputArchive ar(ss);
        ar(*this);

    } else {

        cereal::PortableBinaryOutputArchive ar(ss);
        ar(*this);
        const auto& str = ss.str();
        VLOG(9) << "Wrote " << str.length() << " bytes" << endl;
        R__b.WriteStdString(str);
    }
}

ostream& TEvent::Print(ostream& s) const {
    if(Reconstructed)
        s << "=== Reconstructed:\n" << *Reconstructed;
    if(MCTrue)
        s << "=== MCTrue:\n" << *MCTrue;
    return s;
}
