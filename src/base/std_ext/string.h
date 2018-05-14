#pragma once

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

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

inline bool string_starts_with(std::string const& value, std::string const& beginning)
{
    if (beginning.size() > value.size()) return false;
    return std::equal(beginning.begin(), beginning.end(), value.begin());
}

/**
 * @brief replace a substring
 * @param subject is modified
 * @param search
 * @param replace
 */
inline void replace(std::string& subject,
                    const std::string& search,
                    const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}

inline std::string replace_str(const std::string& input,
                               const std::string& search,
                               const std::string& replace) {
    std::string subject = input;
    std_ext::replace(subject, search, replace);
    return subject;
}


/**
 * @brief string_sanitize removes nasty whitespace from string
 * @param in
 * @return
 */
inline std::string string_sanitize(std::string s) {
  constexpr auto whitespace = " \t\f\v\n\r";
  auto start = s.find_first_not_of(whitespace);
  auto end = s.find_last_not_of(whitespace);
  if(start == std::string::npos && end == std::string::npos)
    return "";
  s.erase(0, start);
  s.erase((end - start) + 1);
  return s;
}

inline bool contains(const std::string& str, const std::string& needle) {
    return str.find(needle) != str.npos;
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
 * @brief Remove all occurences of charachter ch from str
 * @param str String to modify
 * @param ch character to remove
 */
inline void remove_char(std::string& str, const char ch)
{
    str.erase(std::remove(str.begin(), str.end(), ch), str.end());
}

/**
 * @brief Remove all occurences of certain charachters from str
 * @param str String to modify
 * @param chars initializer_list containing characters to be removed
 */
inline void remove_chars(std::string& str, const std::initializer_list<const char> chars)
{
    for (const auto ch : chars)
        remove_char(str, ch);
}

inline std::string basename(const std::string& filenamepath) {

    auto pos = filenamepath.find_last_of("/");

    if(pos==filenamepath.npos) {
        return filenamepath;
    } else {
        return filenamepath.substr(pos,filenamepath.size()-pos);
    }


}

inline std::vector<std::string> tokenize_string(const std::string& str, const std::string& delim) {
    std::vector<std::string> tokens;
    std::string::size_type p = 0;
    std::string::size_type np = 0;

    do {
        np = str.find(delim, p);
        tokens.push_back(str.substr(p, np==str.npos? np : np-p));
        p = np+delim.size();
    } while(np != str.npos);

    return tokens;
}

inline std::string concatenate_string(const std::vector<std::string>& tokens, const std::string& delim) {
    std::stringstream ss;
    for(auto it = tokens.begin(); it != tokens.end(); ++it) {
        ss << *it;
        if(it != std::prev(tokens.end()))
            ss << delim;
    }
    return ss.str();
}

inline std::string to_lower(const std::string& str) {
    std::string lower = str;
    std::transform(str.begin(), str.end(), lower.begin(),
                   [](unsigned char c){ return std::tolower(c); }
    );
    return lower;
}

inline std::string to_upper(const std::string& str) {
    std::string upper = str;
    std::transform(str.begin(), str.end(), upper.begin(),
                   [](unsigned char c){ return std::toupper(c); }
    );
    return upper;
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

    std::string str() const { return stream_.str(); }

    operator std::string () const { return stream_.str(); }

    formatter(const formatter &) = delete;
    formatter & operator = (formatter &) = delete;

private:
    std::stringstream stream_;

};

}} // namespace ant::std_ext
