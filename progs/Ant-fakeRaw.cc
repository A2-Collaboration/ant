#include "expconfig/ExpConfig.h"
#include "analysis/input/ant/AntReader.h"
#include "unpacker/UnpackerAcqu.h"
#include "unpacker/detail/UnpackerAcqu_legacy.h"
#include "unpacker/RawFileReader.h"

#include "expconfig/setups/Setup.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"
#include "tree/TAntHeader.h"

#include "base/WrapTFile.h"
#include "base/std_ext/system.h"
#include "base/CmdLine.h"
#include "base/Logger.h"

#include <memory>
#include <signal.h>
#include <fstream>

using namespace std;
using namespace ant;
using namespace ant::unpacker;

static volatile bool interrupt = false;



int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) {
        interrupt = true;
    });

    TCLAP::CmdLine cmd("Ant-fakeRaw", ' ', "0.1");

    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup manually by name",false,"","setup");
    auto cmd_input  = cmd.add<TCLAP::ValueArg<string>>("i","input","Input files",true,"","filename");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",true,"","filename");
    auto cmd_headerfile = cmd.add<TCLAP::ValueArg<string>>("","header","Acqu Mk2 file to for header extraction",true,"","acqufile");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    // check if header file is readable
    const string& headerfile = cmd_headerfile->getValue();
    {
        string errmsg;
        if(!std_ext::system::testopen(headerfile, errmsg)) {
            LOG(ERROR) << "Cannot open headerfile '" << headerfile << "': " << errmsg;
            return EXIT_FAILURE;
        }
    }

    // check if input file is readable
    const string& inputfile = cmd_input->getValue();
    {
        string errmsg;
        if(!std_ext::system::testopen(inputfile, errmsg)) {
            LOG(ERROR) << "Cannot open inputfile '" << inputfile << "': " << errmsg;
            return EXIT_FAILURE;
        }
    }

    auto inputrootfile = make_shared<WrapTFileInput>(inputfile);

    // check if there's a previous AntHeader present,
    // which could tell us the SetupName
    TAntHeader* previous_AntHeader = nullptr;
    if(inputrootfile->GetObject<TAntHeader>("AntHeader",previous_AntHeader)) {
        const auto& setupname = previous_AntHeader->SetupName;
        if(!setupname.empty()) {
            ExpConfig::Setup::SetManualName(setupname);
            LOG(INFO) << "Setup name set to '" << setupname << "' from input file";
        }
        else
            LOG(WARNING) << "Found AntHeader in input files, but SetupName was empty";
    }

    // override the setup name from cmd line
    if(cmd_setup->isSet()) {
        const auto& setupname = cmd_setup->getValue();
        ExpConfig::Setup::SetManualName(setupname, false);
        LOG(INFO) << "Commandline override setup name to '" << setupname << "'";
    }


    auto setup = ExpConfig::Setup::GetLastFound();
    if(!setup) {
        LOG(ERROR) << "No Setup found. Maybe specify one with --setup?";
        return EXIT_FAILURE;
    }
    // get mapping from setup
    std::vector<UnpackerAcquConfig::hit_mapping_t> hit_mappings;
    {
        std::vector<UnpackerAcquConfig::scaler_mapping_t> scaler_mappings;
        auto config = dynamic_pointer_cast<UnpackerAcquConfig, ExpConfig::Setup>(setup);
        if(!config) {
            LOG(ERROR) << "Provided setup does not know how to unpack Acqu data";
            return EXIT_FAILURE;
        }
        config->BuildMappings(hit_mappings, scaler_mappings);
    }



    const auto& outputfilename = cmd_output->getValue();
    ofstream outputfile(outputfilename);

    // copy the header from given file
    {
        std::vector<uint32_t> buffer;
        RawFileReader r;
        r.open(headerfile);
        if(!r) {
            LOG(ERROR) << "Cannot open headerfile for reading";
            return EXIT_FAILURE;
        }
        // header is usually contained within the first 32k bytes
        r.expand_buffer(buffer, 0x8000/sizeof(uint32_t));
        buffer.resize(10*0x8000/sizeof(uint32_t));
        outputfile.write(reinterpret_cast<const char*>(buffer.data()),
                         buffer.size()*sizeof(uint32_t));
    }

    analysis::input::AntReader reader(inputrootfile, nullptr, nullptr);

    unsigned nEvents = 0;
    TEvent event;
    // ensure first databuffer has Mk2 data marker
    std::vector<uint32_t> databuffer{acqu::EMk2DataBuff};
    while(reader.ReadNextEvent(event)) {
        if(interrupt)
            break;

        vector<uint32_t> eventbuffer;
        {
            for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
                auto it_hit_mapping = std::find_if(hit_mappings.begin(), hit_mappings.end(),
                                                [&readhit] (const UnpackerAcquConfig::hit_mapping_t& m) {
                    return readhit.DetectorType == m.LogicalChannel.DetectorType &&
                           readhit.ChannelType == m.LogicalChannel.ChannelType &&
                           readhit.Channel == m.LogicalChannel.Channel;
                });
                if(it_hit_mapping == hit_mappings.end()) {
                    LOG(WARNING) << "Did not find mapping for given readhit";
                    continue;
                }
                if(it_hit_mapping->RawChannels.size() != 1) {
                    LOG(WARNING) << "Ignoring mapping with non-trivial raw channel";
                    continue;
                }
                auto& rawchannel = it_hit_mapping->RawChannels.front();
                if(rawchannel.NoMask() != rawchannel.Mask) {
                    LOG(WARNING) << "Ignoring readhit with non-trivial mask";
                    continue;
                }
                // handle multihits
                auto nHits = readhit.RawData.size()/sizeof(uint16_t);
                for(unsigned i=0;i<nHits;i++) {
                    auto ptr = reinterpret_cast<const uint16_t*>(readhit.RawData.data())+i;
                    eventbuffer.emplace_back();
                    auto acquhit = reinterpret_cast<acqu::AcquBlock_t*>(addressof(eventbuffer.back()));
                    acquhit->id = rawchannel.RawChannel;
                    acquhit->adc = *ptr;
                }
            }

            // finish with end-event marker
            eventbuffer.push_back(acqu::EEndEvent);
        }

        // each eventbuffer is filled into the buffer with its eventID, eventlength (in bytes)
        // and some possible end-of-databuffer marker, in total 3 extra words maximum
        if(databuffer.size()+eventbuffer.size()+3 > 10*0x8000/sizeof(uint32_t)) {
            databuffer.emplace_back(acqu::EBufferEnd);
            databuffer.resize(10*0x8000/sizeof(uint32_t));
            outputfile.write(reinterpret_cast<const char*>(databuffer.data()),
                             databuffer.size()*sizeof(uint32_t));
            databuffer.clear();
            databuffer.emplace_back(acqu::EMk2DataBuff);
        }

        // fill event into buffer
        databuffer.emplace_back(nEvents);
        databuffer.emplace_back(eventbuffer.size()*sizeof(uint32_t));
        databuffer.insert(databuffer.end(), eventbuffer.begin(), eventbuffer.end());

        nEvents++;
    }

    // dump last databuffer
    databuffer.emplace_back(acqu::EBufferEnd);
    databuffer.resize(10*0x8000/sizeof(uint32_t));
    outputfile.write(reinterpret_cast<const char*>(databuffer.data()),
                     databuffer.size()*sizeof(uint32_t));

    // write EEndBuffer
    databuffer.clear();
    databuffer.emplace_back(acqu::EEndBuff);
    databuffer.emplace_back(acqu::EBufferEnd);
    databuffer.resize(10*0x8000/sizeof(uint32_t));
    outputfile.write(reinterpret_cast<const char*>(databuffer.data()),
                     databuffer.size()*sizeof(uint32_t));

    LOG(INFO) << nEvents << " events processed";

    outputfile.close();

    return EXIT_SUCCESS;
}
