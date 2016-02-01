#pragma once

#include "SlowControl.h"
#include "tree/TSlowControl.h"

#include <list>

namespace ant {
namespace analysis {
namespace input {

/**
 * @brief Fills the values from a TSlowControl into the corresponding field of a SlowControl object
 * @param slc SlowControl object fo fill
 * @param value The TSlowControl object
 */
void FillSlowControl(SlowControl& slc, const TSlowControl& value);

/**
 * @brief Get the list of requested TSlowControl::Keys from a SlowControl object.
 * @param slc
 * @return list of keys
 */
std::list<TSlowControl::Key> RequestedKeys(const SlowControl& slc);

}
}
}
