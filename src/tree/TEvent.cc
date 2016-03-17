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

TEvent::TEvent(const TID& id_reconstructed)
{
    MakeReconstructed(id_reconstructed);
}

TEvent::TEvent(const TID& id_reconstructed, const TID& id_mctrue)
{
    MakeReconstructedMCTrue(id_reconstructed, id_mctrue);
}

void TEvent::MakeReconstructed(const TID& id_reconstructed)
{
    reconstructed = std_ext::make_unique<TEventData>(id_reconstructed);
}

void TEvent::MakeMCTrue(const TID& id_mctrue)
{
    mctrue = std_ext::make_unique<TEventData>(id_mctrue);
}

void TEvent::MakeReconstructedMCTrue(const TID& id_reconstructed, const TID& id_mctrue)
{
    MakeReconstructed(id_reconstructed);
    MakeMCTrue(id_mctrue);
}

void TEvent::ClearDetectorReadHits()
{
    if(reconstructed)
        Reconstructed().ClearDetectorReadHits();
    if(mctrue)
        MCTrue().ClearDetectorReadHits();
}

void TEvent::EnsureTempBranches()
{
    if(!reconstructed) {
        MakeReconstructed(TID());
        empty_reconstructed = true;
    }
    if(!mctrue) {
        MakeMCTrue(TID());
        empty_mctrue = true;
    }
}

void TEvent::ClearTempBranches()
{
    if(empty_reconstructed) {
        reconstructed = nullptr;
        empty_reconstructed = false;
    }
    if(empty_mctrue) {
        mctrue = nullptr;
        empty_mctrue = false;
    }
}

TEvent& TEvent::operator=(TEvent&&) = default;
TEvent::TEvent(TEvent&&) = default;


ostream& TEvent::Print(ostream& s) const
{
    if(reconstructed)
        s << "> Reconstructed:\n" << Reconstructed();
    if(mctrue)
        s << "> MCTrue:\n" << MCTrue();
    return s;
}



