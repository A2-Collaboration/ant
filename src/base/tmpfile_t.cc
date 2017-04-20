#include "tmpfile_t.h"


#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <iterator>
#include <fstream>
#include <sstream>

using namespace std;
using namespace ant;

size_t tmpfile_t::tmpfiles = 0;

tmpfile_t::tmpfile_t() : filename(), testdata()
{
  // obtain some random filename
  char filename_[128];
  strcpy(filename_, "anttmpfile.XXXXXX");
  if(mkstemp(filename_) == -1)
    throw runtime_error("Cannot create tmpfile");
  filename = filename_;
}

tmpfile_t::tmpfile_t(const tmpfolder_t& tmpfolder, const string& extension)
{
    filename = tmpfolder.foldername + "/" + to_string(tmpfiles) + extension;
    tmpfiles++;
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
    // might fail...
    remove(filename.c_str());
}


tmpfolder_t::tmpfolder_t()
{
    // obtain some random foldername
    char foldername_[128];
    strcpy(foldername_, "anttmpfile.XXXXXX");
    if(mkdtemp(foldername_) == NULL)
      throw runtime_error("Cannot create tmpfolder");
    foldername = foldername_;
}

tmpfolder_t::~tmpfolder_t()
{
    stringstream cmd;
    cmd << "rm -r " << foldername << " 2>/dev/null";
    auto unused = system(cmd.str().c_str());
    ++unused; // workaround for unused result / variable
}
