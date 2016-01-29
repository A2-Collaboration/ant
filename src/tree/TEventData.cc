#include "TEventData.h"
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

ostream& TEventData::Print(ostream& s) const {
    return s << "TEventData ID=" << ID;
}

void TEventData::Streamer(TBuffer& R__b) {

    stringstream ss;

    if (R__b.IsReading()) {

        ID.Streamer(R__b);
        DetectorHits.Streamer(R__b);
        Tagger.Streamer(R__b);

        string s;
        R__b.ReadStdString(s);
        VLOG(9) << "Read  " << s.length() << " bytes" << endl;
        ss << s;
        cereal::PortableBinaryInputArchive ar(ss);
        ar(*this);

    } else {

        ID.Streamer(R__b);
        DetectorHits.Streamer(R__b);
        Tagger.Streamer(R__b);


        cereal::PortableBinaryOutputArchive ar(ss);
        ar(*this);
        const auto& str = ss.str();
        VLOG(9) << "Wrote " << str.length() << " bytes" << endl;
        R__b.WriteStdString(str);
    }
}
