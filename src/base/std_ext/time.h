#pragma once

#include "string.h"

#include <ctime>


namespace ant {
namespace std_ext {

inline std::string ctime(const time_t& time) {
  // std::ctime returns some carriage return
  return string_sanitize(std::ctime(std::addressof(time)));
}

inline time_t to_time_t(const std::string& str, bool is_end) {
  // try to parse it as ISO standard format
  std::tm tm_str;
  strptime(str.c_str(), "%Y-%m-%d", &tm_str);
  // set the fields not defined by the given string
  tm_str.tm_isdst = 0;
  tm_str.tm_sec = is_end ? 59 : 0;
  tm_str.tm_min = is_end ? 59 : 0;
  tm_str.tm_hour = is_end ? 23 : 0;
  return mktime(&tm_str);
}

inline bool time_before(const time_t& moment, const std::string& before) {
  return moment<to_time_t(before, true);
}

inline bool time_after(const time_t& moment, const std::string& after) {
  return moment>to_time_t(after, false);
}

inline bool time_between(const time_t& moment,
                         const std::string& begin,
                         const std::string& end) {
  return time_before(moment, end) && time_after(moment, begin);
}

}} // namespace ant::std_ext
