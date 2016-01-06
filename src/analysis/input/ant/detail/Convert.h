#pragma once

#include "analysis/data/Candidate.h"

#include <memory>

namespace ant {

struct TEvent;
struct TCandidate;
struct TTaggerHit;
struct TCluster;

namespace analysis {

namespace data {
    struct Cluster;
    struct Event;
    struct TaggerHit;
}

namespace input {

/**
 * @brief Convert a TEvent from the reconstruct stage into an analysis style event
 * @param event
 * @return a shared ptr to the new event
 */

struct Converter {

    static data::Event Convert(const TEvent& event);

    static data::CandidatePtr Convert(const TCandidate& candidate);

    static data::TaggerHit Convert(const TTaggerHit& taggerhit);

    static data::Cluster Convert(const TCluster& cluster);

};


}
}
}
