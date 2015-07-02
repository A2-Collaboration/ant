

// this ensure the correct forward 
// declaration of lzma_stream
#include <lzma.h>
#define RAWFILEREADER_H_IMPL
#include "RawFileReader.h"
#undef RAWFILEREADER_H_IMPL

#include <cstdio> // for BUFSIZ

using namespace std;
using namespace ant;




RawFileReader::XZ::XZ(const std::string &filename) : 
  Plain(filename),
  decompressFailed(false),
  gcount_(0)
{
  strm = (lzma_stream*)malloc(sizeof(lzma_stream));
  try {  
    init_decoder();
  }
  catch(Exception e) {
    cleanup();
    // propagate the exception
    throw e;
  }
  strm->next_in = NULL;
  strm->avail_in = 0;  
}

RawFileReader::XZ::~XZ() {
  cleanup();
}

void RawFileReader::XZ::cleanup()
{
  lzma_end(strm);
  free(strm);
}

void RawFileReader::XZ::read(char* s, streamsize n) { 
 
  constexpr streamsize inbufsiz = BUFSIZ;  
    
  lzma_action action = LZMA_RUN;
  
  char inbuf[inbufsiz];
  
  strm->next_out = reinterpret_cast<uint8_t*>(s);
  strm->avail_out = n;
  
  while (true) {
    
    if (strm->avail_in == 0 && !Plain::eof()) {
      strm->next_in = reinterpret_cast<uint8_t*>(inbuf);
      Plain::read(inbuf, inbufsiz); // read from the underlying file into buffer
      strm->avail_in = Plain::gcount();
      
      if (!Plain::operator bool()) {
        throw Exception("Error while reading from compressed input file");
      }
      
      if (Plain::eof())
        action = LZMA_FINISH;
    }
    
    lzma_ret ret = lzma_code(strm, action);
    
    if (ret == LZMA_STREAM_END) {
      gcount_ = n - strm->avail_out;
      decompressFailed = true; // incomplete decompress
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





