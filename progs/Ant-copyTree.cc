
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/WrapTFile.h"

#include "TTree.h"

#include <iostream>

using namespace std;
using namespace ant;

int main (int argc, char** argv)
{

    SetupLogger();

    TCLAP::CmdLine cmd("Ant-copyTree - Copy all trees from input to output, up to max events", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","filename");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",true,"","filename");
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    WrapTFileInput inputfile(cmd_input->getValue());

    auto intrees = inputfile.GetListOf<TTree>();

    if(intrees.empty()) {
        LOG(ERROR) << "No trees found in inputfile";
        return EXIT_FAILURE;
    }

    WrapTFileOutput outputfile(cmd_output->getValue());
    outputfile.cd();

    long long maxevents = cmd_maxevents->isSet()
            ? cmd_maxevents->getValue().back()
            :  numeric_limits<long long>::max();

    for(auto in : intrees) {
        VLOG(1) << "Found tree " << in->GetName() << " with entries=" << in->GetEntries();
        in->CloneTree(maxevents);
    }

    return EXIT_SUCCESS;
}
