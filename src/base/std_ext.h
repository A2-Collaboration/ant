#ifndef ANT_STD_EXT_H
#define ANT_STD_EXT_H

#include <string>
#include <sstream>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <ctime>
#include <iostream>

#include <cxxabi.h> // works only with GCC?

namespace ant {

namespace std_ext {

inline void insertRange(std::vector<unsigned>& v, unsigned start, unsigned stop) {
  int length = stop-start+1;
  if(length<1)
    return;
  std::vector<unsigned> v_(static_cast<size_t>(length));
  std::iota(v_.begin(), v_.end(), start);
  v.insert(v.end(), v_.cbegin(), v_.cend());
}

template<typename T>
inline bool contains(const std::vector<T> v, const T& val) {
  return std::find(v.cbegin(), v.cend(), val) != v.cend();
}

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

template<typename T>
std::string getTypeAsString() {
  return abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);
}

inline double degree_to_radian(const double degree) {
  return degree * M_PI / 180.0;
}

// copied from C++14 draft
// https://isocpp.org/files/papers/N3656.txt
// once we use C++14 by default, we can remove this here
template<class T> struct _Unique_if {
  typedef std::unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]> {
  typedef std::unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]> {
  typedef void _Known_bound;
};

template<class T, class... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n) {
  typedef typename std::remove_extent<T>::type U;
  return std::unique_ptr<T>(new U[n]());
}

template<class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args&&...) = delete;


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

} // namespace ant::std_ext

} // namespace ant



#endif // ANT_STD_EXT
