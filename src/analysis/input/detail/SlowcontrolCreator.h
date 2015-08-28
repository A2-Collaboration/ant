#pragma once

#include "data/Slowcontrol.h"
#include "tree/TSlowControl.h"

#include <list>

namespace ant {
namespace analysis {

/**
 * @brief Fills the values from a TSlowControl into the corresponding field of a Slowcontrol object
 * @param slc Slowcontrol object fo fill
 * @param value The TSlowcontrol object
 */
void FillSlowcontrol(data::Slowcontrol& slc, const TSlowControl& value);

/**
 * @brief Get the list of requested TSlowControl::Keys from a Slowcontrol object.
 * @param slc
 * @return list of keys
 */
std::list<TSlowControl::Key> RequestedKeys(const data::Slowcontrol& slc);

}
}
