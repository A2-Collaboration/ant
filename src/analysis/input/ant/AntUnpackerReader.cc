#include "AntUnpackerReader.h"

#include "data/Event.h"
#include "detail/Convert.h"
#include "base/ReadTFiles.h"
#include "base/std_ext.h"
#include "tree/TEvent.h"

#include "TTree.h"

#include <memory>
#include <stdexcept>

using namespace std;
using namespace ant;
using namespace ant::input;

AntUnpackerReader::AntUnpackerReader(
        unique_ptr<Unpacker::Reader> unpacker_reader,
        unique_ptr<tree::UnpackerWriter> unpacker_writer) :
    reader(move(unpacker_reader)),
    writer(move(unpacker_writer))
{

}



AntUnpackerReader::~AntUnpackerReader() {}

std::shared_ptr<Event> AntUnpackerReader::ReadNextEvent()
{
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
