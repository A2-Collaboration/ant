#include "expconfig/ExpConfig.h"
#include "unpacker/Unpacker.h"
#include "unpacker/RawFileReader.h"

#include "reconstruct/Reconstruct.h"

#include "tree/TDetectorRead.h"
#include "expconfig/setups/SetupRegistry.h"

#include "base/WrapTFile.h"
#include "base/std_ext/system.h"
#include "base/CmdLine.h"
#include "base/Logger.h"

#include "TTree.h"

#include <memory>
#include <signal.h>

using namespace std;
using namespace ant;

static volatile bool keepRunning = true;

int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) {
        keepRunning = false;
    });

    TCLAP::CmdLine cmd("Ant-rawdump", ' ', "0.1");

    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_input  = cmd.add<TCLAP::ValueArg<string>>("i","input","Input files",true,"","filename");
    auto cmd_ADCs = cmd.add<TCLAP::MultiArg<string>>("a","adc","Ranges of Acqu ADC numbers, e.g. 400-412;5;10",true,"AcquADCs");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    if(std_ext::system::isInteractive())
        RawFileReader::OutputPerformanceStats = 3;

    // check if input file is readable
    const string& inputfile = cmd_input->getValue();
    string errmsg;
    if(!std_ext::system::testopen(inputfile, errmsg)) {
        LOG(ERROR) << "Cannot open inputfile '" << inputfile << "': " << errmsg;
        return 1;
    }


    // construct the ADC option for Setup_Raw
    std::shared_ptr<OptionsList> setup_opts = make_shared<OptionsList>();
    stringstream AcquADC_option;
    AcquADC_option << "AcquADC=";
    for(auto adcstr : cmd_ADCs->getValue()) {
        AcquADC_option << adcstr << " ";
    }
    setup_opts->SetOption(AcquADC_option.str());
    ant::expconfig::SetupRegistry::SetSetupOptions(setup_opts);


    // override the setup name
    ExpConfig::Setup::ManualName = "Setup_Raw";

    // now we can try to open the files with an unpacker
    std::unique_ptr<Unpacker::Reader> unpacker = nullptr;

    try {
        unpacker = Unpacker::Get(inputfile);
    }
    catch(Unpacker::Exception e) {
        LOG(ERROR) << "Unpacker exception: " << e.what();
        return 1;
    }
    catch(RawFileReader::Exception e) {
        LOG(ERROR) << "Unpacker: Error opening file "<<inputfile<<": " << e.what();
        return 1;
    }
    catch(...) {
        LOG(ERROR) << "Cannot create unpacker for input file";
        return 1;
    }

    // the real output file, create it here to get all
    // further ROOT objects into this output file
    unique_ptr<WrapTFileOutput> masterFile;
    unique_ptr<vector<unsigned>> ADC_index;
    unique_ptr<vector<unsigned>> ADC_value;
    TTree* tree = nullptr;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                           WrapTFileOutput::mode_t::recreate,
                                                           true); // cd into masterFile upon creation

        tree = masterFile->CreateInside<TTree>("AcquRawADC","AcquRawADC");
        ADC_index = std_ext::make_unique<vector<unsigned>>();
        ADC_value = std_ext::make_unique<vector<unsigned>>();
        tree->Branch("ADC_index", "std::vector<unsigned int>", ADC_index.get());
        tree->Branch("ADC_value", "std::vector<unsigned int>", ADC_value.get());
    }

    // run the loop
    while(auto item = unpacker->NextItem()) {
        if(!keepRunning)
            break;
        TDetectorRead* detread = dynamic_cast<TDetectorRead*>(item.get());
        if(detread==nullptr)
            continue;

        if(!masterFile) {
            LOG(INFO) << *detread;
            continue;
        }

        ADC_index->resize(0);
        ADC_value->resize(0);
        for(const TDetectorReadHit& readhit : detread->Hits) {
            // assume raw data to be 16bit chunks
            const auto& rawData = readhit.RawData;
            if(rawData.size() % 2 != 0)
                continue;
            for(size_t i=0;i<rawData.size()/2;i++) {
                const uint16_t* rawVal = reinterpret_cast<const uint16_t*>(addressof(rawData[2*i]));
                ADC_index->push_back(readhit.Channel);
                ADC_value->push_back(*rawVal);
            }
        }

        tree->Fill();
    }

    return 0;
}
