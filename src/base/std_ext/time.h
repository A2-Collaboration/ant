#pragma once

#include "base/std_ext/string.h"

#include <cstring>
#include <ctime>


namespace ant {
namespace std_ext {

inline std::string to_iso8601(const time_t& time) {
    char buf[sizeof "2011-10-08T07:07:09Z"];
    strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", gmtime(&time));
    // std::ctime returns some carriage return
    return buf;
}

inline std::tm to_tm(const std::string& str, const std::string& fmt) {
    std::tm tm_str;
    strptime(str.c_str(), fmt.c_str(), std::addressof(tm_str));
    return tm_str;
}

inline time_t to_time_t(const std::tm& tm_) {
    tm tm_copy(tm_);
    return std::mktime(std::addressof(tm_copy));
}

inline time_t to_time_t(const std::string& str, bool is_end) {
  // try to parse it as ISO standard format
  std::tm tm_str = to_tm(str, "%Y-%m-%d");
  // set the fields not defined by the given string
  tm_str.tm_sec = is_end ? 59 : 0;
  tm_str.tm_min = is_end ? 59 : 0;
  tm_str.tm_hour = is_end ? 23 : 0;
  return std::mktime(std::addressof(tm_str));
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

inline int get_day_of_week(int day, int month, int year) {
    // Use Zeller's formula with convention from tm struct,
    // return days since Sunday
    month += 1;
    year += 1900;
    if (month < 3) {
       year--;
       month += 12;
    }
    int dayno = 1 + day + 2*month + 3*(month + 1)/5 + year + year/4 - year/100 + year/400;
    dayno %= 7;
    return dayno;
}

inline bool guess_dst(std::tm& time) {
    // try to use mktime to guess it
    time.tm_isdst = -1;
    std::mktime(std::addressof(time));
    // but catch this one hour in october
    // because it's implementation specific what happens then
    if(time.tm_mon != 9 || time.tm_hour != 2)
        return true;
    // get last sunday in october
    int lastsunday = 0;
    for(int day=1;day<=31;day++) {
        if(get_day_of_week(day, time.tm_mon, time.tm_year) == 0)
            lastsunday = day;
    }
    return time.tm_mday != lastsunday;
}

}} // namespace ant::std_ext
