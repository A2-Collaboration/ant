#include "Logger.h"

#include "TError.h"
#include <sstream>

// setup the logger, will be compiled as a little library
INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace ant::logger;

void SetupLogger(int argc, char* argv[]) {
    START_EASYLOGGINGPP(argc, argv);
    SetupLogger();
}

void SetupLogger() {
    el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
    el::Configurations loggerConf;
    loggerConf.setToDefault();
    loggerConf.setGlobally(el::ConfigurationType::Format, "%datetime [%level] %fbase:%line : %msg");
    loggerConf.setGlobally(el::ConfigurationType::ToFile, "false");
    loggerConf.set(el::Level::Verbose,  el::ConfigurationType::Format, "%datetime [%level-%vlevel] %fbase:%line : %msg");
    el::Loggers::reconfigureLogger("default", loggerConf);

    // set ROOT error handler
    SetErrorHandler([] (
                    int level, Bool_t, const char *location,
                    const char *msg) {
        // Bool_t is abort, not used for now...

        if (level < gErrorIgnoreLevel)
            return;

        // somehow ROOT issues a very strange warning
        // when using gROOT->FindObjectAny (called in ant::canvas::FindTCanvas to search for TCanvas)
        if(level == kWarning &&
           string(location) == "TClass::TClass" &&
           string(msg) == "no dictionary for class iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&> is available")
            return; // ignore it for now

        if(level == kWarning &&
           string(location) == "TTree::Bronch" &&
           string(msg) == "Using split mode on a class: TLorentzVector with a custom Streamer")
            return; // ignore it for now

        if(level == kWarning &&
           string(location) == "TTree::Bronch" &&
           string(msg) == "Using split mode on a class: TVector2 with a custom Streamer")
            return; // ignore it for now


        stringstream ss;
        ss << "ROOT<" << location << ">: " << msg;
        if(level < kInfo) { LOG(INFO) << ss.str(); }
        else if(level == kInfo) { LOG(INFO) << ss.str(); }
        else if(level == kWarning) { LOG(WARNING) << ss.str(); }
        else {
            LOG(ERROR) << ss.str();
            if(DebugInfo::nProcessedEvents>=0)
                LOG(DEBUG) << "nProcessEvents=" << DebugInfo::nProcessedEvents;
            if(DebugInfo::nUnpackedBuffers>=0)
                LOG(DEBUG) << "nUnpackedBuffers=" << DebugInfo::nUnpackedBuffers;
        }

    });
}

long long DebugInfo::nProcessedEvents = -1;
int DebugInfo::nUnpackedBuffers = -1;

