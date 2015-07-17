#ifndef CONVERT_H
#define CONVERT_H

#include <memory>

class ant::Event;
class TEvent;

namespace ant {
namespace input {
namespace ant {

/**
 * @brief Convert a TEvent from the reconstruct stage into an analysis style event
 * @param event
 * @return a shared ptr to the new event
 */
std::shared_ptr<ant::Event> Convert(const std::shared_ptr<TEvent>& event);

}
}
}

#endif
