
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "Unpacker.h"

#include "Logger.h"

using namespace std;
using namespace ant;

int main(int argc, char* argv[]) {
  START_EASYLOGGINGPP(argc, argv);
  el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, "%datetime [%level] %fbase : %msg");

  auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_7892.dat.xz");

  unpacker->NextItem();

  //cout << unpacker->OpenFile("scratch/CBTaggTAPS_7892.dat.xz") << endl;

//  cout << unpacker->NextItem() << endl;

//  vector<int> v{1,2,3,4};

//  LOG(INFO) << v;

//  for(size_t i=0;i<100;i++) {
//    LOG_N_TIMES(10,WARNING) << "OHOH";
//  }

  return EXIT_SUCCESS;
}
