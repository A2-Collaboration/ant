#include "tmpfile_t.h"


#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <iterator>
#include <fstream>

using namespace std;

tmpfile_t::tmpfile_t() : filename(), testdata()
{
  // obtain some random filename
  char filename_[128];
  strcpy(filename_, "anttmpfile.XXXXXX");
  if(mkstemp(filename_) == -1)
    throw runtime_error("Cannot create tmpfile "+filename);
  filename = filename_;
}


void tmpfile_t::write_testdata() {
  ofstream outfile(filename);
  ostream_iterator<uint8_t> outiterator(outfile);
  copy(testdata.begin(), testdata.end(), outiterator);
  outfile.close();
  if(!outfile)
    throw runtime_error("Cannot write data to "+filename);
}


tmpfile_t::~tmpfile_t() {
  // do not really
  remove(filename.c_str());
}
