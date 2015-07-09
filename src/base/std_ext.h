#ifndef ANT_STD_EXT_H
#define ANT_STD_EXT_H

#include <string>
#include <sstream>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <ctime>

namespace ant {

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

inline std::string ctime(const time_t& time) {
  // std::ctime returns some carriage return
  return string_sanitize(std::ctime(std::addressof(time)));
}


template<typename From, typename To>
inline std::unique_ptr<To> static_cast_uptr(std::unique_ptr<From>&& ptr) {
  static_assert(std::is_base_of<To, From>::value, "Type From does not derive from type To");
  return std::unique_ptr<To>(static_cast<To*>(ptr.release()));
}

/**
 * @brief The formatter class
 *
 * Used to create formatted strings with stringstream
 */
class formatter
{
public:
    formatter() {}
    ~formatter() {}

    template <typename Type>
    formatter & operator << (const Type & value)
    {
        stream_ << value;
        return *this;
    }

    std::string str() const         { return stream_.str(); }
    operator std::string () const   { return stream_.str(); }

private:
    std::stringstream stream_;

    formatter(const formatter &);
    formatter & operator = (formatter &);
};

} // namespace ant::std_ext

} // namespace ant



#endif // ANT_STD_EXT
