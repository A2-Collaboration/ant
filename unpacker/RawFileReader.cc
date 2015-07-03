

// this ensure the correct forward
// declaration of lzma_stream
#include <lzma.h>
#define RAWFILEREADER_H_IMPL
#include "RawFileReader.h"
#undef RAWFILEREADER_H_IMPL

#include <cstdio> // for BUFSIZ
#include <cstring> // for strerror

using namespace std;
using namespace ant;


RawFileReader::XZ::XZ(const std::string &filename, const size_t inbufsize) :
  PlainBase(filename),
  inbuf(inbufsize),
  decompressFailed(false),
  gcount_(0),
  eof_(false)
{
  strm = (lzma_stream*)malloc(sizeof(lzma_stream));
  if(strm==nullptr)
    throw Exception("Cannot allocate memory");
  try {
    init_decoder();
  }
  catch(Exception e) {
    cleanup();
    // propagate the exception
    throw e;
  }
}

RawFileReader::XZ::~XZ() {
  cleanup();
}

void RawFileReader::XZ::cleanup()
{
  lzma_end(strm);
  free(strm);
}

void RawFileReader::XZ::init_decoder()
{
  *strm = LZMA_STREAM_INIT;

  lzma_ret ret = lzma_stream_decoder(
        strm, UINT64_MAX, LZMA_CONCATENATED);

  // Return successfully if the initialization went fine.
  if (ret == LZMA_OK) {
    return;
  }

  // otherwise throw some exceptions
  switch (ret) {
  case LZMA_MEM_ERROR:
    throw Exception("Memory allocation failed");
  case LZMA_OPTIONS_ERROR:
    throw Exception("Unsupported decompressor flags");
  default:
    throw Exception("Unknown error, possibly a bug");
  }
}

void RawFileReader::XZ::read(char* s, streamsize n) {

  lzma_action action = PlainBase::eof() ? LZMA_FINISH : LZMA_RUN;

  strm->next_out = reinterpret_cast<uint8_t*>(s);
  strm->avail_out = n;

  while (true) {

    if (strm->avail_in == 0 && !PlainBase::eof()) {
      // read from the underlying compressed file into buffer
      PlainBase::read(reinterpret_cast<char*>(inbuf.data()), inbuf.size());
      strm->next_in = inbuf.data();
      strm->avail_in = PlainBase::gcount();

      if(PlainBase::eof()) {
        action = LZMA_FINISH;
      }
      else if (!PlainBase::operator bool()) {
        throw Exception(string("Error while reading from compressed input file: ")
                        +string(strerror(errno)));
      }
    }

    lzma_ret ret = lzma_code(strm, action);

    if(ret == LZMA_STREAM_END) {
      gcount_ = n - strm->avail_out; // number of decompressed bytes
      if(strm->avail_out > 0)
        eof_ = true;
      return;
    }

    if(strm->avail_out == 0) {
      gcount_ = n;
      return;
    }

    if(ret == LZMA_OK)
      continue;

    decompressFailed = true;

    // ret indicates some error condition
    switch (ret) {
    case LZMA_MEM_ERROR:
      throw Exception("Memory allocation failed");
      break;

    case LZMA_FORMAT_ERROR:
      throw Exception("The input is not in the .xz format");
      break;

    case LZMA_OPTIONS_ERROR:
      throw Exception("Unsupported compression options");
      break;

    case LZMA_DATA_ERROR:
      throw Exception("Compressed file is corrupt");
      break;

    case LZMA_BUF_ERROR:
      throw Exception("Compressed file is truncated or "
                      "otherwise corrupt");

    default:
      throw Exception("Unknown error, possibly a bug");

    }
  }
}








