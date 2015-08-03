#pragma once

#include <memory>

namespace ant {
class Cluster;
class Event;
class Candidate;
class TaggerHit;

class TEvent;
class TCandidate;
class TTaggerHit;
class TCluster;

namespace input {

/**
 * @brief Convert a TEvent from the reconstruct stage into an analysis style event
 * @param event
 * @return a shared ptr to the new event
 */

struct Converter {

static Event Convert(const TEvent& event);

static std::shared_ptr<ant::Candidate> Convert(const TCandidate& candidate);

static std::shared_ptr<ant::TaggerHit> Convert(const TTaggerHit& taggerhit);

static Cluster Convert(const TCluster& cluster);

};


}
}

