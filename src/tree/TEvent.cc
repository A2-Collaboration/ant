#include "TEvent.h"
#include "TEventData.h"
#include "stream_TBuffer.h"

#include "base/cereal/types/polymorphic.hpp"
#include "base/cereal/types/memory.hpp"
#include "base/cereal/types/vector.hpp"
#include "base/cereal/types/list.hpp"                 // hidden in TParticleTree_t ...
#include "base/cereal/archives/binary.hpp"


#include "base/std_ext/memory.h"
#include "base/Logger.h"

#include "TClass.h"

#include <streambuf>

using namespace std;
using namespace ant;

// use some versioning
CEREAL_CLASS_VERSION(TEvent, ANT_TEVENT_VERSION)

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

// create some TBuffer to std::streambuf interface





void TEvent::Streamer(TBuffer& R__b)
{
    stream_TBuffer buf(R__b);
    iostream inoutstream(addressof(buf));

    if (R__b.IsReading()) {
        cereal::BinaryInputArchive ar(inoutstream);
        ar(*this);
    }
    else {
        cereal::BinaryOutputArchive ar(inoutstream);
        ar(*this);
    }
}

ostream& TEvent::Print(ostream& s) const {
    if(Reconstructed)
        s << "> Reconstructed:\n" << *Reconstructed;
    if(MCTrue)
        s << "> MCTrue:\n" << *MCTrue;
    return s;
}

std::unique_ptr<TEvent> TEvent::MakeReconstructed(const TID& id) {
    auto event = std_ext::make_unique<TEvent>();
    event->Reconstructed = std_ext::make_unique<TEventData>(id);
    return event;
}

void TEvent::ClearDetectorReadHits()
{
    if(Reconstructed)
        Reconstructed->ClearDetectorReadHits();
    if(MCTrue)
        MCTrue->ClearDetectorReadHits();
}

TEvent::TEvent() : Reconstructed(), MCTrue() {}
TEvent::~TEvent() {}


