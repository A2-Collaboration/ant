#include "AntReader.h"

#include "data/Event.h"

#include "detail/Convert.h"

#include "tree/UnpackerWriter.h"
#include "tree/TEvent.h"
#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"

#include "base/std_ext.h"
#include "base/Logger.h"

#include "TTree.h"

#include <memory>
#include <stdexcept>

using namespace std;
using namespace ant;
using namespace ant::input;

AntReader::AntReader(
        unique_ptr<Unpacker::Reader> unpacker_reader,
        std::unique_ptr<Reconstruct_traits> reconstruct
        ) :
    reader(move(unpacker_reader)),
    writer(nullptr),
    reconstruct(move(reconstruct)),
    haveReconstruct(false),
    nEvents(0),
    writeUncalibrated(false),
    writeCalibrated(false)
{

}

AntReader::~AntReader() {}

void AntReader::EnableUnpackerWriter(
        const string& outputfile,
        bool uncalibratedDetectorReads,
        bool calibratedDetectorReads
        )
{
    writer = std_ext::make_unique<tree::UnpackerWriter>(outputfile);
    writeUncalibrated = uncalibratedDetectorReads;
    writeCalibrated = calibratedDetectorReads;
    LOG(INFO) << "Writing unpacker stage output to " << outputfile;
    if(writeUncalibrated)
        LOG(INFO) << "Write UNcalibrated detectors reads (BEFORE DoReconstruct) to " << outputfile;
    if(writeCalibrated)
        LOG(INFO) << "Write calibrated detectors (AFTER DoReconstruct) reads to " << outputfile;
}

bool AntReader::ReadNextEvent(Event& event, TSlowControl&)
{
    while(auto item = reader->NextItem()) {
        // we use ROOT's machinery to identify derived class types,
        // because it's much faster than dynamic_cast (but also potentially unsafe)
        const TClass* isA = item->IsA();

        if(isA == THeaderInfo::Class()) {
            const THeaderInfo* headerInfo = reinterpret_cast<THeaderInfo*>(item.get());
            if(reconstruct) {
                reconstruct->Initialize(*headerInfo);
                haveReconstruct = true;
                LOG(INFO) << "Found THeaderInfo in unpacker datastream, initialized Reconstruct";
            }
        }
        else if(isA == TDetectorRead::Class()) {
            TDetectorRead* detread = reinterpret_cast<TDetectorRead*>(item.get());

            if(writer && writeUncalibrated)
                writer->Fill(item);

            if(haveReconstruct) {
                auto tevent = reconstruct->DoReconstruct(*detread);
                event = input::Convert(*tevent);
                if(writer) {
                    if(writeCalibrated)
                        writer->Fill(item);
                    writer->Fill(move(tevent));
                }
                return true;
            }

            nEvents++;
            // skip the writing of the detector read item
            // because item is handled above
            continue;
        }
        else if(isA == TEvent::Class()) {
            const TEvent* tevent = reinterpret_cast<TEvent*>(item.get());
            event = input::Convert(*tevent);
            return true;
        }

        // by default, we write the items to the file
        if(writer)
            writer->Fill(item);
    }

    return false;
}

