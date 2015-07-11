
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "Unpacker.h"
#include "expconfig/ExpConfig.h"
#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"
#include "unpacker/RawFileReader.h"

#include "base/Logger.h"
#include "base/Format.h"

#include <chrono>

using namespace std;
using namespace ant;

int main(int argc, char* argv[]) {
  SetupLogger(argc, argv);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();

  //std::vector<double> v{2,3};
  //string str = fmt::format("Bla {}", v);
  //cout << str << endl;

  //auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_9227.dat");
  //auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_7892.dat");
  auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_5711.dat.xz");
//  auto unpacker = Unpacker::Get("scratch/oneevent-small.dat");
  unsigned nReads = 0;

  while(auto item = unpacker->NextItem()) {
    //cout << *item << endl;
    const TDetectorRead* read = dynamic_cast<const TDetectorRead*>(item.get());
    if(read == nullptr)
      continue;
    //cout << read << endl;
    nReads++;
    if(nReads % 10000 == 0) {
      VLOG(5) << "Unpacked " << nReads << " detector reads";
    }
  }

  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  cout << "Analyzed " << nReads/elapsed_seconds.count() << " Reads/s" << endl;

//  THeaderInfo header(TDataRecord::ID_t(0,0), 0, "", 0);
//  auto config = shared_ptr<ExpConfig::Module>(ExpConfig::Get(header));
//  LOG(INFO) << "Got config: " << config << endl;

  //cout << unpacker->OpenFile("scratch/CBTaggTAPS_7892.dat.xz") << endl;

//  cout << unpacker->NextItem() << endl;

//  vector<int> v{1,2,3,4};

//  LOG(INFO) << v;

//  for(size_t i=0;i<100;i++) {
//    LOG_N_TIMES(10,WARNING) << "OHOH";
//  }

  return EXIT_SUCCESS;
}
