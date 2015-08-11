#include "tmpfile_t.h"


#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <iterator>
#include <fstream>
#include <cassert>

using namespace std;
using namespace ant;

tmpfile_t::tmpfile_t(const string& extension) : filename(), testdata()
{
    string templ = "anttmpfile.XXXXXX" + extension;
    assert(templ.size() <128);

    // obtain some random filename
    char filename_[128];
    strcpy(filename_, templ.c_str());
    if(mkstemps(filename_,extension.size()) == -1)
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
