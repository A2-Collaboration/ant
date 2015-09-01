#pragma once

#include <string>
#include <sstream>
#include <algorithm>

namespace ant {
namespace std_ext {

/**
 * @brief Check if a string ends with a substr. For example useful for file extensions.
 * @param value String to check
 * @param ending Ending to check for
 * @return true it value ends with ending
 */
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

/**
 * @brief Remove all occurences of substr from str
 * @param str String to modify
 * @param substr substring to remove
 */
inline void removesubstr(std::string& str, const std::string& substr) {

    std::string::size_type pos = 0;

    while(true) {

        pos = str.find(substr, pos);

        if(pos == str.npos)
            break;

        str.erase(pos, substr.length());
    }
}

/**
 * @brief The formatter class
 *
 * Used to create formatted strings with stringstream
 */
class formatter
{
public:
    formatter() : stream_() {}
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

}} // namespace ant::std_ext
