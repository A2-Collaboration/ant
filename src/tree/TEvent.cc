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


// teach cereal some ROOT types

template<class Archive>
void save(Archive& archive,
          const TVector3& m)
{
    archive( m.X(), m.Y(), m.Z() );
}

template<class Archive>
void load(Archive & archive,
          TVector3& m)
{
    double x,y,z;
    archive(x, y, z);
    m.SetXYZ(x,y,z);
}

template<class Archive>
void save(Archive& archive,
          const TLorentzVector& m)
{
    archive( m.Px(), m.Py(), m.Pz(), m.E() );
}

template<class Archive>
void load(Archive & archive,
          TLorentzVector& m)
{
    double px,py,pz,e;
    archive(px, py, pz, e);
    m.SetPxPyPzE(px,py,pz,e);
}

// tell cereal to use the correct TParticle load/save due to inheritance from TLorentzVector

namespace cereal
{
  template <class Archive>
  struct specialize<Archive, TParticle, cereal::specialization::member_load_save> {};
}


void TEvent::Streamer(TBuffer& R__b) {

    /// \todo create streambuf class from TBuffer to improve performance?
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
        s << "> Reconstructed:\n" << *Reconstructed;
    if(MCTrue)
        s << "> MCTrue:\n" << *MCTrue;
    return s;
}

ostream& TEvent::Data::Print(ostream& s) const {
    s << "ID=" << ID << endl;

    s << ">> DetectorReadHits" << endl;
    for(auto& i: DetectorReadHits)
        s << i << endl;

    s << ">> SlowControls" << endl;
    for(auto& i : SlowControls)
        s << i << endl;

    s << ">> UnpackerMessages" << endl;
    for(auto& i : UnpackerMessages)
        s << i << endl;

    s << ">> Tagger" << endl << Tagger;

    s << ">> Clusters" << endl;
    for(auto& i : Clusters)
        s << *i << endl;

    s << ">> Candidates" << endl;
    for(auto& i : Candidates)
        s << *i << endl;

    s << ">> Particles" << endl;
    for(auto& i : Particles)
        s << *i << endl;

    return s;
}
