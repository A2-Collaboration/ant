#include "Logger.h"
// setup the logger, will be compiled as a little library
INITIALIZE_EASYLOGGINGPP

void SetupLogger(int argc, char* argv[]) {
  START_EASYLOGGINGPP(argc, argv);
  el::Configurations loggerConf;
  loggerConf.setToDefault();
  loggerConf.setGlobally(el::ConfigurationType::Format, "%datetime [%level] %fbase:%line : %msg");
  loggerConf.set(el::Level::Verbose,  el::ConfigurationType::Format, "%datetime [%level-%vlevel] %fbase:%line : %msg");
  el::Loggers::reconfigureLogger("default", loggerConf);
}