#pragma once

#include <memory>

namespace ant {
class TEvent;
class TCandidate;
class TTaggerHit;

class Event;
class Track;
class TaggerHit;

namespace input {

/**
 * @brief Convert a TEvent from the reconstruct stage into an analysis style event
 * @param event
 * @return a shared ptr to the new event
 */
std::shared_ptr<ant::Event> Convert(const TEvent& event);

std::shared_ptr<ant::Track> Convert(const TCandidate& track);

std::shared_ptr<ant::TaggerHit> Convert(const TTaggerHit& taggerhit);
}
}

