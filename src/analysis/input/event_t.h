#pragma once

#include "tree/TEvent.h"

namespace ant {
namespace analysis {
namespace input {

struct event_t : TEvent {

    // we're a wrapper around TEvent
    event_t() = default;
    explicit event_t(TEvent event) :
        TEvent(std::move(event))
    {}

    bool empty_reconstructed = false;
    bool empty_mctrue = false;

    bool HasReconstructed() const { return reconstructed!=nullptr; }
    bool HasMCTrue() const { return mctrue!=nullptr; }

    void MakeReconstructed(const TID& id_reconstructed);
    void MakeMCTrue(const TID& id_mctrue);
    void MakeReconstructedMCTrue(const TID& id_reconstructed, const TID& id_mctrue);

    void ClearDetectorReadHits();
    void EnsureTempBranches();
    void ClearTempBranches();
};

}}}