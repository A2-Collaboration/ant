
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
  el::Loggers::reconfigureLogger("default",
                                 el::ConfigurationType::Format,
                                 "%datetime [%level] %fbase : %msg");

  el::Configurations loggerConf;
  loggerConf.setToDefault();
  loggerConf.setGlobally(el::ConfigurationType::Format, "%datetime [%level] %fbase : %msg");
  loggerConf.set(el::Level::Verbose,  el::ConfigurationType::Format, "%datetime [%level-%vlevel] %fbase : %msg");
  el::Loggers::reconfigureLogger("default", loggerConf);

  auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_7892.dat.xz");

  LOG(INFO) << "Got item: " << unpacker->NextItem();

  //cout << unpacker->OpenFile("scratch/CBTaggTAPS_7892.dat.xz") << endl;

//  cout << unpacker->NextItem() << endl;

//  vector<int> v{1,2,3,4};

//  LOG(INFO) << v;

//  for(size_t i=0;i<100;i++) {
//    LOG_N_TIMES(10,WARNING) << "OHOH";
//  }

  return EXIT_SUCCESS;
}
