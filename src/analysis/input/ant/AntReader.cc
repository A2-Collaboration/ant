#include "AntReader.h"

#include "data/Event.h"

#include "detail/Convert.h"

#include "tree/UnpackerWriter.h"
#include "tree/TEvent.h"
#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"
#include "tree/TSlowControl.h"

#include "base/Logger.h"

#include "TTree.h"

#include <memory>
#include <stdexcept>

using namespace std;
using namespace ant;
using namespace ant::analysis::input;
using namespace ant::analysis::data;

AntReader::AntReader(
        unique_ptr<Unpacker::Reader> unpacker_reader,
        std::unique_ptr<Reconstruct_traits> reconstruct
        ) :
    reader(move(unpacker_reader)),
    reconstruct(move(reconstruct)),
    haveReconstruct(false),
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
    if(writeCalibrated && writeUncalibrated)
        throw Exception("Writing calibrated AND uncalibrated detector reads in one file makes no sense");

    if(writeUncalibrated)
        LOG(INFO) << "Write UNcalibrated detectors reads (BEFORE DoReconstruct) to " << outputfile;
    else if(writeCalibrated)
        LOG(INFO) << "Write calibrated detectors (AFTER DoReconstruct) reads to " << outputfile;
    else
        LOG(INFO) << "Writing unpacker TEvents to " << outputfile;
}

unique_ptr<TSlowControl> AntReader::ReadNextSlowControl()
{
    return move(buffered_slowcontrol);
}

bool AntReader::ReadNextEvent(Event& event)
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
                writer->Fill(item.get());

            if(haveReconstruct) {
                MemoryPool<TEvent>::Item tevent = reconstruct->DoReconstruct(*detread);
                event = Converter::Convert(*tevent);
                if(writer) {
                    if(writeCalibrated)
                        writer->Fill(item.get());
                    else if(!writeUncalibrated)
                        writer->Fill(tevent.get());
                }
                return true;
            }

            // skip the writing of the detector read item
            // because writing is already handled above
            continue;
        }
        else if(isA == TEvent::Class()) {
            const TEvent* tevent = reinterpret_cast<TEvent*>(item.get());
            event = Converter::Convert(*tevent);
            // do not write TEvent's again to the file
            return true;
        }
        else if(isA == TSlowControl::Class()) {
            if(writer)
                writer->Fill(item.get());
            buffered_slowcontrol = unique_ptr<TSlowControl>(reinterpret_cast<TSlowControl*>(item.release()));
            return false;
        }
        /// @todo handle TUnpackerMessage, especially DAQ errors might be of interest for the analysis

        // by default, we write the items to the file
        if(writer)
            writer->Fill(item.get());
    }

    return false;
}

