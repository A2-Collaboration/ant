#include "AntReader.h"

#include "tree/TEvent.h"

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
    writeUncalibrated = uncalibratedDetectorReads;
    writeCalibrated = calibratedDetectorReads;
    if(writeCalibrated && writeUncalibrated)
        throw Exception("Writing calibrated AND uncalibrated detector reads in one file makes no sense");

    if(writeUncalibrated)
        LOG(INFO) << "Write UNcalibrated detectors reads (BEFORE DoReconstruct) to " << outputfile;
    else if(writeCalibrated)
        LOG(INFO) << "Write calibrated detectors (AFTER DoReconstruct) reads to " << outputfile;
    else
        LOG(INFO) << "Writing unpacker Events to " << outputfile;
}

double AntReader::PercentDone() const
{
    return reader->PercentDone();
}

bool AntReader::ReadNextEvent(TEvent& event)
{
    auto eventptr = reader->NextEvent();
    /// \todo Handle reading/writing here?

    if(eventptr) {
        event = move(*eventptr);
        return true;
    }

    return false;

//    while() {
//        // we use ROOT's machinery to identify derived class types,
//        // because it's much faster than dynamic_cast (but also potentially unsafe)
//        const TClass* isA = item->IsA();

//        if(isA == TDetectorRead::Class()) {
//            TDetectorRead* detread = reinterpret_cast<TDetectorRead*>(item.get());

//            if(writer && writeUncalibrated)
//                writer->Fill(item.get());

//            if(reconstruct) {
//                MemoryPool<TEvent>::Item tevent = reconstruct->DoReconstruct(*detread);
//                event = Converter::Convert(*tevent);
//                if(writer) {
//                    if(writeCalibrated)
//                        writer->Fill(item.get());
//                    else if(!writeUncalibrated)
//                        writer->Fill(tevent.get());
//                }
//                return true;
//            }

//            // skip the writing of the detector read item
//            // because writing is already handled above
//            continue;
//        }
//        else if(isA == Event::Class()) {
//            const Event* tevent = reinterpret_cast<Event*>(item.get());
//            event = Converter::Convert(*tevent);
//            // do not write Event's again to the file
//            return true;
//        }
//        else if(isA == TSlowControl::Class()) {
//            if(writer)
//                writer->Fill(item.get());
//            buffered_slowcontrol = unique_ptr<TSlowControl>(reinterpret_cast<TSlowControl*>(item.release()));
//            return false;
//        }
//        /// @todo handle TUnpackerMessage, especially DAQ errors might be of interest for the analysis

//        // by default, we write the items to the file
//        if(writer)
//            writer->Fill(item.get());
//    }

//    return false;
}

