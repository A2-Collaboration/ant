#include "catch.hpp"
#include "RawFileReader.h"
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <iterator>
#include <cstring>
#include <cstdio>
#include <exception>

using namespace std;

// some configs
constexpr streamsize totalSize = 10*BUFSIZ;
constexpr streamsize chunkSize = totalSize/17;
constexpr streamsize inbufSize = BUFSIZ;

// use some little class which cleans up after itself properly
struct tmpfile_t {
  string filename;
  vector<uint8_t> testdata;

  tmpfile_t() {
    // obtain some random filename
    char filename_[128];
    strcpy(filename_, "rawfilereader.XXXXXX");
    if(mkstemp(filename_) == -1)
      throw runtime_error("Cannot create tmpfile for test");
    filename = filename_;
  }


  void write_testdata() {
    ofstream outfile(filename);
    ostream_iterator<uint8_t> outiterator(outfile);
    copy(testdata.begin(), testdata.end(), outiterator);
    outfile.close();
    REQUIRE(outfile); // file properly written and closed
  }

  ~tmpfile_t() {
    if(remove(filename.c_str()) != 0)
       throw runtime_error("Cannot cleanup tmpfile for test");
  }
};

void dotest(bool, streamsize, streamsize, streamsize);
void doendianness();

TEST_CASE("Test RawFileReader: nocompress, one chunk", "[unpacker]") {
  dotest(false, totalSize, totalSize, inbufSize);
}

TEST_CASE("Test RawFileReader: compress, one chunk", "[unpacker]") {
  dotest(true, totalSize, totalSize, inbufSize);
}

TEST_CASE("Test RawFileReader: nocompress, chunks", "[unpacker]") {
  dotest(false, totalSize, chunkSize, inbufSize);
}

TEST_CASE("Test RawFileReader: compress, chunks", "[unpacker]") {
  dotest(true, totalSize, chunkSize, inbufSize);
}

TEST_CASE("Test RawFileReader: weird stuff", "[unpacker]") {
  dotest(true, 100, 7, 40); // inputbuffer smaller than output buffers
}

TEST_CASE("Test RawFileReader: uint32_t endianness","[unpacker]") {
  doendianness();
}

void doendianness() {
  tmpfile_t f;

  f.testdata = {0x10, 0x20, 0x30, 0x40};
  f.write_testdata();
  ant::RawFileReader reader;
  REQUIRE_NOTHROW(reader.open(f.filename));

  vector<uint32_t> indata(1);
  REQUIRE_NOTHROW(reader.read(indata.data(), indata.size()));

  for(size_t i=0;i<f.testdata.size();i++) {
    const uint32_t byte = f.testdata[i];
    const uint32_t word = (indata[0] >> 8*i) & 0xff;
    REQUIRE(byte == word);
  }
}

void dotest(bool compress,
            streamsize totalSize,
            streamsize chunkSize,
            streamsize inbufSize) {
  tmpfile_t f;
  // write some testdata to given temporary filename
  f.testdata.resize(totalSize);
  //int n(0);
  //generate(testdata.begin(), testdata.end(), [&n] {return n++;});
  generate(f.testdata.begin(), f.testdata.end(), rand);
  f.write_testdata();

  // make a little detour for compression
  // the RawFileReader should be able to decompress it
  // transparently
  if(compress) {
    //compress it first
    const string& xz_cmd = string("xz ")+f.filename;
    REQUIRE(system(xz_cmd.c_str()) == 0);
    f.filename += ".xz"; // xz changes the filename
  }

  // then continue reading in the file with the RawFileReader
  ant::RawFileReader reader;
  REQUIRE_NOTHROW(reader.open(f.filename, inbufSize));

  { // INFO scope
    INFO("Check if file is open");
    REQUIRE(reader); // uses the operator bool() of RawFileReader
  }

  // read in the data again and compare bytewise
  vector<uint8_t> indata(f.testdata.size());

  // two different reading modes are possible
  if(totalSize == chunkSize) {
    REQUIRE_NOTHROW(reader.read((char*)&indata[0], indata.size()));
    REQUIRE(reader.gcount()==indata.size());
    REQUIRE(!reader.eof());
    // reading just one more byte should not be possible
    // and should be indicated by eof
    REQUIRE_NOTHROW(reader.read((char*)&indata[0], 1));
    REQUIRE(reader.gcount()==0);
    REQUIRE(reader.eof());
  }
  else {
    const size_t nChunks = totalSize / chunkSize;
    size_t offset = 0;
    for(size_t i=0;i<nChunks;i++) {
      REQUIRE_NOTHROW(reader.read((char*)&indata[offset], chunkSize));
      REQUIRE(reader.gcount()==chunkSize);
      REQUIRE(!reader.eof());
      offset += chunkSize;
    }
    // read one additional (possibly incomplete) chunk
    REQUIRE_NOTHROW(reader.read((char*)&indata[offset], chunkSize));
    REQUIRE(reader.gcount() == totalSize-nChunks*chunkSize);
    REQUIRE(reader.gcount() <= chunkSize); // never more than requested

    // the eof flag is there when trying to read past the end of a file
    // so, if by conincidence, the chunks fitted into totalSize with remainder,
    // the r.eof is not set yet
    if(reader.gcount() == chunkSize)
      REQUIRE(!reader.eof());
    else
      REQUIRE(reader.eof());

    // reading just one more byte should not be possible
    // and should be indicated by eof
    REQUIRE_NOTHROW(reader.read((char*)&indata[0], 1));
    REQUIRE(reader.gcount()==0);
    REQUIRE(reader.eof());
  }

  // test the read in data
  const bool inputEqualsOutput = indata == f.testdata;
  REQUIRE(inputEqualsOutput);
}



