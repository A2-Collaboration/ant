#pragma once

#include "base/bitflag.h"

namespace ant {
namespace analysis {
namespace input {

enum class reader_flag_t {
    IsSource,                  // reader act as TEvent source (does not only amend them)
    ProvidesSlowControl,     // reader provides slow control information
};

using reader_flags_t = bitflag<reader_flag_t>;

}}}