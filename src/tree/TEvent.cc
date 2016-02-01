#include "TEvent.h"
#include "TClass.h"

#include "base/cereal/types/polymorphic.hpp"
#include "base/cereal/types/memory.hpp"
#include "base/cereal/types/vector.hpp"
#include "base/cereal/types/list.hpp"                 // hidden in TParticleTree_t ...
#include "base/cereal/archives/portable_binary.hpp"


#include "base/std_ext/memory.h"
#include "base/Logger.h"

#include <sstream>
#include <streambuf>
#include <iostream>

using namespace std;
using namespace ant;

const TParticleList TEvent::Data::PTypeList::empty;

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

/// \todo the implementation below causes a segfault at weird places
/// when writing larger amounts of data,
/// probably TBuffer is not treated correctly...

//class stream_TBuffer : public std::streambuf {
//public:
//    explicit stream_TBuffer(TBuffer& tbuffer_) :
//        tbuffer(tbuffer_)
//    {
//        if(tbuffer.IsReading()) {
//            // reading uses the default std::streambuf behaviour
//            const auto begin = tbuffer.Buffer()+tbuffer.Length();
//            const auto end = tbuffer.Buffer()+tbuffer.BufferSize();
//            setg(begin, begin, end);
//        }
//    }
//private:
//    // the streambuf interface for writing
//    streamsize xsputn(const char_type* s, streamsize n) override {
//        tbuffer.WriteFastArray(s, n);
//        return n;
//    }

//    int_type overflow(int_type ch) override {
//        if(ch != traits_type::eof()) {
//            tbuffer.WriteChar(ch);
//        }
//        return ch;
//    }

//    // forbid copy
//    stream_TBuffer(const stream_TBuffer&) = delete;
//    stream_TBuffer& operator=(const stream_TBuffer&) = delete;

//    // hold a reference to the buffer for writing business
//    TBuffer& tbuffer;
//};



void TEvent::Streamer(TBuffer& R__b)
{
//    stream_TBuffer buf(R__b);
//    iostream inoutstream(addressof(buf));

    stringstream ss;
    if (R__b.IsReading()) {
        string s;
        R__b.ReadStdString(s);
        ss << s;
        cereal::PortableBinaryInputArchive ar(ss);
        ar(*this);
    }
    else {
        cereal::PortableBinaryOutputArchive ar(ss);
        ar(*this);
        R__b.WriteStdString(ss.str());
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
    event->Reconstructed = std_ext::make_unique<TEvent::Data>(id);
    return event;
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
    for(auto& i : Particles.GetAll())
        s << *i << endl;

    return s;
}
