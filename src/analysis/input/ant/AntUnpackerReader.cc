#include "AntUnpackerReader.h"

#include "data/Event.h"

#include "detail/Convert.h"


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

AntUnpackerReader::AntUnpackerReader(
        unique_ptr<Unpacker::Reader> unpacker_reader,
        std::unique_ptr<Reconstruct_traits> reconstruct
        ) :
    reader(move(unpacker_reader)),
    writer(nullptr),
    reconstruct(move(reconstruct)),
    nEvents(0),
    writeUncalibrated(false),
    writeCalibrated(false)
{

}



AntUnpackerReader::~AntUnpackerReader() {}

std::shared_ptr<Event> AntUnpackerReader::ReadNextEvent()
{
    while(auto item = reader->NextItem()) {


        // we use ROOT's machinery to identify derived class types,
        // because it's much faster than dynamic_cast (but also potentially unsafe)
        const TClass* isA = item->IsA();

        if(isA == TDetectorRead::Class()) {
            TDetectorRead* detread = reinterpret_cast<TDetectorRead*>(item.get());

            if(writer && writeUncalibrated)
                writer->Fill(item);

            if(reconstruct) {
                auto event = reconstruct->DoReconstruct(*detread);
                if(writer) {
                    if(writeCalibrated)
                        writer->Fill(item);
                    writer->Fill(event);
                }
                return input::Convert(*event);
            }

            nEvents++;
            // skip the writing of the detector read item
            // because item is handled above
            continue;
        }
        else if(isA == TEvent::Class()) {
            const TEvent* event = reinterpret_cast<TEvent*>(item.get());
            return input::Convert(*event);
        }

        // by default, we write the items to the file
        if(writer)
            writer->Fill(item);
    }

    return nullptr;
}

bool AntUnpackerReader::hasData() const {
    return false;
}

long long AntUnpackerReader::EventsRead() const
{
    return 0;
}

long long AntUnpackerReader::TotalEvents() const
{
    return -1;
}
