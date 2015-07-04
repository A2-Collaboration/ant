

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

void RawFileReader::open(const string &filename, const size_t inbufsize) {
  // open it as plain raw file
  ifstream file(filename.c_str());

  // already check here if the open was ok
  if(!file)
    throw Exception(string("Error when opening file ")
                    +filename
                    +": "
                    +string(strerror(errno)));

  if(XZ::test(file)) {
    p = std::unique_ptr<PlainBase>(new XZ(filename, inbufsize));
  }
  else {
    p = std::unique_ptr<PlainBase>(new PlainBase(filename));
  }
}


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

bool RawFileReader::XZ::test(ifstream& file) {
  // check some magic bytes in the beginning
  // to determine possible compression
  std::vector<char> file_bytes(6);
  file.seekg(0, file.beg); // ensure start of file
  file.read(file_bytes.data(), file_bytes.size());


  vector<char> magic_bytes_xz{ static_cast<char>(0xFD), '7', 'z', 'X', 'Z', 0x00 };

  return file_bytes == magic_bytes_xz;
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
