#pragma once

#include <memory>
#include "analysis/data/Cluster.h"

namespace ant {
class TEvent;
class TCandidate;
class TTaggerHit;

class Event;
class Candidate;
class TaggerHit;
class TCluster;

namespace input {

/**
 * @brief Convert a TEvent from the reconstruct stage into an analysis style event
 * @param event
 * @return a shared ptr to the new event
 */
ant::Event Convert(const TEvent& event);

std::shared_ptr<ant::Candidate> Convert(const TCandidate& candidate);

std::shared_ptr<ant::TaggerHit> Convert(const TTaggerHit& taggerhit);

ant::Cluster Convert(const TCluster& cluster);
}
}

