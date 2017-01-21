#include "event_t.h"

#include "tree/TEventData.h"
#include "base/std_ext/memory.h"

using namespace ant;
using namespace ant::analysis::input;

void event_t::MakeReconstructed(const TID& id_reconstructed)
{
    reconstructed = std_ext::make_unique<TEventData>(id_reconstructed);
}

void event_t::MakeMCTrue(const TID& id_mctrue)
{
    mctrue = std_ext::make_unique<TEventData>(id_mctrue);
}

void event_t::MakeReconstructedMCTrue(const TID& id_reconstructed, const TID& id_mctrue)
{
    MakeReconstructed(id_reconstructed);
    MakeMCTrue(id_mctrue);
}

void event_t::ClearDetectorReadHits()
{
    if(HasReconstructed())
        Reconstructed().ClearDetectorReadHits();
    if(HasMCTrue())
        MCTrue().ClearDetectorReadHits();
}

void event_t::EnsureTempBranches()
{
    if(!HasReconstructed()) {
        MakeReconstructed(TID());
        empty_reconstructed = true;
    }
    if(!HasMCTrue()) {
        MakeMCTrue(TID());
        empty_mctrue = true;
    }
}

void event_t::ClearTempBranches()
{
    if(empty_reconstructed) {
        reconstructed = nullptr;
        empty_reconstructed = false;
    }
    if(empty_mctrue) {
        mctrue = nullptr;
        empty_mctrue = false;
    }
}
