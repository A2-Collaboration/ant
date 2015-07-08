#ifndef FORMAT_H
#define FORMAT_H

#include "detail/format.h"
#include <vector>

namespace fmt {

// really ugly way of implementing it
template<typename T>
std::string format_vector(StringRef format_str, const std::vector<T> v) {
  switch(v.size()) {
  case 0:
    return format_str;
  case 1:
    return format(format_str, v[0]);
  case 2:
    return format(format_str, v[0], v[1]);
  case 3:
    return format(format_str, v[0], v[1], v[2]);
  case 4:
    return format(format_str, v[0], v[1], v[2], v[3]);
  default:
    throw FormatError("Too many vector elements to format it");
  }
}

}

#endif // FORMAT_H
