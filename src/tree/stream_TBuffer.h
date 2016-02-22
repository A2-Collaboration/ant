#pragma once

#include "TBuffer.h"
#include <streambuf>

namespace ant {

class stream_TBuffer : public std::streambuf {
public:
    explicit stream_TBuffer(TBuffer& tbuffer_) :
        tbuffer(tbuffer_)
    {
        if(tbuffer.IsReading()) {
            // reading uses the default std::streambuf behaviour
            const auto begin = tbuffer.Buffer()+tbuffer.Length();
            const auto end = tbuffer.Buffer()+tbuffer.BufferSize();
            setg(begin, begin, end);
        }
    }
private:
    // the streambuf interface for writing
    std::streamsize xsputn(const char_type* s, std::streamsize n) override {
        tbuffer.WriteFastArray(s, n);
        return n;
    }

    int_type overflow(int_type ch) override {
        if(ch != traits_type::eof()) {
            tbuffer.WriteChar(ch);
        }
        return ch;
    }

    // forbid copy
    stream_TBuffer(const stream_TBuffer&) = delete;
    stream_TBuffer& operator=(const stream_TBuffer&) = delete;

    // hold a reference to the buffer for writing business
    TBuffer& tbuffer;
};

}