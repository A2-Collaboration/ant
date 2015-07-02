
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
  RawFileReader r("scratch/CBTaggTAPS_7892.dat.xz");
  if(r) {
    
    vector<uint32_t> buffer(128);
    try {
      r.read((uint8_t*)&buffer[0], 4*buffer.size());
    }
    catch(ant::RawFileReader::Exception e) {
      cerr << "exception " << e.what() << endl;
    }

    if(r) {
      cout << hex;
      for(size_t i=0;i<buffer.size();i++) {
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
