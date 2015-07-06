#include "Logger.h"
// setup the logger, will be compiled as a little library
INITIALIZE_EASYLOGGINGPP

void SetupLogger(int argc, char* argv[]) {
  START_EASYLOGGINGPP(argc, argv);
  el::Loggers::reconfigureLogger("default",
                                 el::ConfigurationType::Format,
                                 "%datetime [%level] %fbase : %msg");

  el::Configurations loggerConf;
  loggerConf.setToDefault();
  loggerConf.setGlobally(el::ConfigurationType::Format, "%datetime [%level] %fbase : %msg");
  loggerConf.set(el::Level::Verbose,  el::ConfigurationType::Format, "%datetime [%level-%vlevel] %fbase : %msg");
  el::Loggers::reconfigureLogger("default", loggerConf);
}