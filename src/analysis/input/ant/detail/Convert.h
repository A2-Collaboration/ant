#ifndef CONVERT_H
#define CONVERT_H

#include <memory>

namespace ant {
class TEvent;
class TTrack;

class Event;
class Track;

namespace input {

/**
 * @brief Convert a TEvent from the reconstruct stage into an analysis style event
 * @param event
 * @return a shared ptr to the new event
 */
std::shared_ptr<ant::Event> Convert(const TEvent &event);

std::shared_ptr<ant::Track> Convert(const TTrack &track);
}
}

#endif
