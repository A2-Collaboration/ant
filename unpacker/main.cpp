
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include "RawFileReader.h"

using namespace std;
using namespace ant;

// fix qtcreator highlighting...
typedef std::uint32_t uint32_t;


int main() {
  RawFileReader r("scratch/CBTaggTAPS_7892.dat");
  if(r) {

    vector<uint32_t> buffer(100);
    unsigned offset = 4;
    unsigned total = 32;

    try {
      r.read((char*)&buffer[0],  4*offset);
      r.read((char*)&buffer[offset],  4*(total-offset));
      //r.read((char*)&buffer[0],  4*8);
    }
    catch(ant::RawFileReader::Exception e) {
      cerr << "exception " << e.what() << endl;
    }

    if(r) {
      cout << hex;
      for(size_t i=0;i<total;i++) {
        if(i>0 && i%8 == 0)
          cout << endl;
        cout << setw(8) << setfill('0') << buffer[i] << " ";
      }
      cout << dec << endl;
    }
    else
      cerr << "Some reading problem" << endl;

  }
  else
    cerr << "File not opened" << endl;
  return EXIT_SUCCESS;
}
