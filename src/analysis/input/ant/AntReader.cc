#include "AntReader.h"

#include "tree/TEvent.h"

#include "base/Logger.h"

#include "TTree.h"

#include <memory>
#include <stdexcept>

using namespace std;
using namespace ant;
using namespace ant::analysis::input;

AntReader::AntReader(
        const std::shared_ptr<WrapTFileInput>& rootfiles,
        unique_ptr<Unpacker::Module> unpacker,
        std::unique_ptr<Reconstruct_traits> reconstruct
        ) :
    unpacker(move(unpacker)),
    reconstruct(move(reconstruct))
{
    /// \todo Implement reading TEvent's from TTree
}

AntReader::~AntReader() {}

double AntReader::PercentDone() const
{
    return unpacker->PercentDone();
}

bool AntReader::ReadNextEvent(TEvent& event)
{
    auto eventptr = unpacker->NextEvent();

    if(eventptr) {
        event = move(*eventptr);
        return true;
    }

    return false;
}

