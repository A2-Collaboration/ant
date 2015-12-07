#pragma once

#include <memory>

namespace ant {

struct TEvent;
struct TCandidate;
struct TTaggerHit;
struct TCluster;

namespace analysis {

namespace data {
    class Cluster;
    class Event;
    class Candidate;
    class TaggerHit;
}

namespace input {

/**
 * @brief Convert a TEvent from the reconstruct stage into an analysis style event
 * @param event
 * @return a shared ptr to the new event
 */

struct Converter {

    static data::Event Convert(const TEvent& event);

    static std::shared_ptr<data::Candidate> Convert(const TCandidate& candidate);

    static data::TaggerHit Convert(const TTaggerHit& taggerhit);

    static data::Cluster Convert(const TCluster& cluster);

};


}
}
}
