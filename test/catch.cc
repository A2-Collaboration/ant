#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "base/Logger.h"


int main( int argc, char* argv[] )
{
  el::Configurations loggerConf;
  loggerConf.setToDefault();
  loggerConf.setGlobally(el::ConfigurationType::Enabled, "false");
  el::Loggers::reconfigureAllLoggers(loggerConf);

  return Catch::Session().run( argc, argv );
}