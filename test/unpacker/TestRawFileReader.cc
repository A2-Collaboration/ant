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

void dotest(bool, streamsize, streamsize);

constexpr streamsize totalSize = 10;
constexpr streamsize chunkSize = 4;

TEST_CASE("Test RawFileReader: nocompress, one chunk", "[unpacker]") {
  dotest(false, totalSize, totalSize);
}

TEST_CASE("Test RawFileReader: compress, one chunk", "[unpacker]") {
  dotest(true, totalSize, totalSize);
}

TEST_CASE("Test RawFileReader: nocompress, chunks", "[unpacker]") {
  dotest(false, totalSize, chunkSize);
}

TEST_CASE("Test RawFileReader: compress, chunks", "[unpacker]") {
  dotest(true, totalSize, chunkSize);
}

void dotest(bool compress, streamsize totalSize, streamsize chunkSize) {

  // use some little class which cleans up after itself properly
  struct tmpfile_t {
    string filename;
    tmpfile_t() {
      // obtain some random filename
      char filename_[128];
      strcpy(filename_, "rawfilereader.XXXXXX");
      if(mkstemp(filename_) == -1)
        throw runtime_error("Cannot create tmpfile for test");
      filename = filename_;
    }

    ~tmpfile_t() {
      if(remove(filename.c_str()) != 0)
         throw runtime_error("Cannot cleanup tmpfile for test");
    }
  };
  tmpfile_t f;

  // write some testdata to given temporary filename
  vector<uint8_t> testdata(totalSize);
  int n(0);
  //generate(testdata.begin(), testdata.end(), [&n] {return n++;});
  generate(testdata.begin(), testdata.end(), rand);
  ofstream outfile(f.filename);
  ostream_iterator<uint8_t> outiterator(outfile);
  copy(testdata.begin(), testdata.end(), outiterator);
  outfile.close();
  REQUIRE(outfile); // file properly written and closed

  // make a little detour for compression
  // the RawFileReader should be able to decompress it
  // transparently
  if(compress) {
    //compress it first
    string xz_cmd = string("xz ")+f.filename;
    REQUIRE(system(xz_cmd.c_str()) == 0);
    f.filename += ".xz"; // xz changes the filename
  }

  // then continue reading in the file with the RawFileReader
  REQUIRE_NOTHROW(ant::RawFileReader r(f.filename));

  ant::RawFileReader r(f.filename);
  { // INFO scope
    INFO("Check if file is open");
    REQUIRE(r); // uses the operator bool() of RawFileReader
  }

  // read in the data again and compare bytewise
  vector<uint8_t> indata(testdata.size());

  // two different reading modes are possible
  if(totalSize == chunkSize) {
    REQUIRE_NOTHROW(r.read((char*)&indata[0], indata.size()));
    REQUIRE(r.gcount()==indata.size());
    REQUIRE(!r.eof());
    // reading just one more byte should not be possible
    // and should be indicated by eof
    REQUIRE_NOTHROW(r.read((char*)&indata[0], 1));
    REQUIRE(r.gcount()==0);
    REQUIRE(r.eof());
  }
  else {
    const size_t nChunks = totalSize / chunkSize;
    size_t offset = 0;
    for(size_t i=0;i<nChunks;i++) {
      REQUIRE_NOTHROW(r.read((char*)&indata[offset], chunkSize));
      REQUIRE(r.gcount()==chunkSize);
      REQUIRE(!r.eof());
      offset += chunkSize;
    }
    // read one additional (possibly incomplete) chunk
    REQUIRE_NOTHROW(r.read((char*)&indata[offset], chunkSize));
    REQUIRE(r.gcount() == totalSize-nChunks*chunkSize);
    REQUIRE(r.gcount() <= chunkSize); // never more than requested

    // the eof flag is there when trying to read past the end of a file
    // so, if by conincidence, the chunks fitted into totalSize with remainder,
    // the r.eof is not set yet
    if(r.gcount() == chunkSize)
      REQUIRE(!r.eof());
    else
      REQUIRE(r.eof());

    // reading just one more byte should not be possible
    // and should be indicated by eof
    REQUIRE_NOTHROW(r.read((char*)&indata[0], 1));
    REQUIRE(r.gcount()==0);
    REQUIRE(r.eof());
  }

  // test the read in data
  const bool inputEqualsOutput = indata == testdata;
  REQUIRE(inputEqualsOutput);
}



