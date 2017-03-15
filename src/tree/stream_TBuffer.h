#pragma once

// ignore warnings from library
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "cereal/cereal.hpp"
#include "cereal/types/polymorphic.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/list.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/bitset.hpp"
#pragma GCC diagnostic pop

#include "TBuffer.h"
#include <streambuf>

namespace ant {

class stream_TBuffer : public std::streambuf {
public:
    explicit stream_TBuffer(TBuffer& tbuffer) :
        tbuffer_(tbuffer)
    {
        if(tbuffer_.IsReading()) {
            // reading uses the default std::streambuf behaviour
            const auto begin = tbuffer_.Buffer()+tbuffer_.Length();
            const auto end = tbuffer_.Buffer()+tbuffer_.BufferSize();
            setg(begin, begin, end);
        }
    }

    // little helper function to call the binary archiver
    // on some class
    template<class T>
    static void DoBinary(TBuffer& tbuffer, T& theClass) {
        stream_TBuffer buf(tbuffer);
        std::iostream inoutstream(addressof(buf));

        if (tbuffer.IsReading()) {
            cereal::BinaryInputArchive ar(inoutstream);
            ar(theClass);
        }
        else {
            cereal::BinaryOutputArchive ar(inoutstream);
            ar(theClass);
        }
    }

private:
    // the streambuf interface for writing
    std::streamsize xsputn(const char_type* s, std::streamsize n) override {
        tbuffer_.WriteFastArray(s, n);
        return n;
    }

    int_type overflow(int_type ch) override {
        if(ch != traits_type::eof()) {
            tbuffer_.WriteChar(ch);
        }
        return ch;
    }

    // forbid copy
    stream_TBuffer(const stream_TBuffer&) = delete;
    stream_TBuffer& operator=(const stream_TBuffer&) = delete;

    // hold a reference to the buffer for writing business
    TBuffer& tbuffer_;
};

}