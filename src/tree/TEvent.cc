#include "TEvent.h"
#include "TEventData.h"
#include "stream_TBuffer.h"

#include "base/std_ext/memory.h"
#include "base/Logger.h"

#include "TClass.h"

#include <streambuf>

using namespace std;
using namespace ant;

// use some versioning
CEREAL_CLASS_VERSION(TEvent, ANT_TEVENT_VERSION)

// tell cereal to use the correct TParticle load/save due to inheritance from LorentzVec

namespace cereal
{
  template <class Archive>
  struct specialize<Archive, TParticle, cereal::specialization::member_load_save> {};
}

// create some TBuffer to std::streambuf interface
void TEvent::Streamer(TBuffer& R__b)
{
    stream_TBuffer::DoBinary(R__b, *this);
}


// other stuff

TEvent::TEvent() : reconstructed(), mctrue() {}
TEvent::~TEvent() {}

TEvent::TEvent(TEvent&&) = default;
TEvent& TEvent::operator=(TEvent&&) = default;


TEvent::TEvent(const TID& id_reconstructed)
{
    reconstructed = std_ext::make_unique<TEventData>(id_reconstructed);
}

TEvent::TEvent(const TID& id_reconstructed, const TID& id_mctrue)
{
    reconstructed = std_ext::make_unique<TEventData>(id_reconstructed);
    mctrue = std_ext::make_unique<TEventData>(id_mctrue);
}

ostream& TEvent::Print(ostream& s) const
{
    if(reconstructed)
        s << "> Reconstructed:\n" << Reconstructed();
    if(mctrue)
        s << "> MCTrue:\n" << MCTrue();
    return s;
}



