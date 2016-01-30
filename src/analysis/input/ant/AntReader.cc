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
    // we expect Reconstructed branch to be filled always
    auto eventptr = unpacker->NextEvent();

    if(eventptr) {
        if(reconstruct)
            reconstruct->DoReconstruct(eventptr->Reconstructed);

        event.Reconstructed = move(eventptr->Reconstructed);

        // the A2Geant unpacker also fills some MCTrue information
        if(eventptr->MCTrue)
            event.MCTrue = move(eventptr->MCTrue);

        return true;
    }

    return false;
}

