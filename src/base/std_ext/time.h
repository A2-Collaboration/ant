#pragma once

#include "string.h"

#include <ctime>


namespace ant {
namespace std_ext {

inline std::string ctime(const time_t& time) {
  // std::ctime returns some carriage return
  return string_sanitize(std::ctime(std::addressof(time)));
}

inline std::tm to_tm(const std::string& str, const std::string& fmt, bool fixdst = true) {
    std::tm tm_str;
    strptime(str.c_str(), fmt.c_str(), std::addressof(tm_str));
    if(fixdst)
        tm_str.tm_isdst = 0;
    return tm_str;
}

inline time_t to_time_t(const std::string& str, bool is_end) {
  // try to parse it as ISO standard format
  std::tm tm_str = to_tm(str, "%Y-%m-%d");
  // set the fields not defined by the given string
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

inline bool is_mest2met_transition(const std::tm& time) {
    // a sunday in october, between 2:00:00 and 2:59:59
    if(time.tm_wday == 0 && time.tm_mon == 9
       && time.tm_hour==2) {
        // determine last sunday of month october
        int lastsunday = 0;
        for(int day=1;day<=31;day++) {
            if(get_day_of_week(day, time.tm_mon, time.tm_year) == 0)
                lastsunday = day;
        }
        if(time.tm_mday == lastsunday)
            return true;
    }
    return false;
}

}} // namespace ant::std_ext
