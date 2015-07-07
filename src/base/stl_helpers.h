#ifndef STL_HELPERS_H
#define STL_HELPERS_H

#include <string>
#include <algorithm>
#include <memory>
#include <type_traits>

namespace std_ext {

inline bool string_ends_with(std::string const& value, std::string const& ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline std::string string_sanitize(const char* in) {
  std::string s(in);
  const std::string& whitespace = " \t\f\v\n\r";
  int start = s.find_first_not_of(whitespace);
  int end = s.find_last_not_of(whitespace);
  if(start == -1 && end == -1)
    return std::string("");
  s.erase(0,start);
  s.erase((end - start) + 1);
  s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());
  return s;
}

template<typename From, typename To>
inline std::unique_ptr<To> static_cast_uptr(std::unique_ptr<From>&& ptr) {
  static_assert(std::is_base_of<To, From>::value, "Type From does not derive from type To");
  return std::unique_ptr<To>(static_cast<To*>(ptr.release()));
}

}

#endif // STL_HELPERS_H
