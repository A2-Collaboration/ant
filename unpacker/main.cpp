
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "Unpacker.h"

using namespace std;
using namespace ant;

// fix qtcreator highlighting...
typedef std::uint32_t uint32_t;


int main() {


  auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_7892.dat.xz");

  cout << unpacker->OpenFile("scratch/CBTaggTAPS_7892.dat.xz") << endl;

  cout << unpacker->NextItem() << endl;

  return EXIT_SUCCESS;
}
