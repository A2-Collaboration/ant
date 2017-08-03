#pragma once

#include <algorithm>
#include <functional>
#include <cmath>
#include <cxxabi.h>
#include <type_traits>
#include <memory>

namespace ant {
namespace std_ext {

using namespace std;

template<class T, class U, class Better>
bool copy_if_better(T& dest, const U& source,
                    Better better) {
    if(std::isfinite(dest) && better(dest, source))
        return false;
    dest = source;
    return true;
}

template<class T, class U>
bool copy_if_greater(T& dest, const U& source) {
    return copy_if_better(dest, source, std::greater<U>());
}


template<class TContainer>
bool begins_with(const TContainer& input, const TContainer& match)
{
    return input.size() >= match.size()
        && equal(match.begin(), match.end(), input.begin());
}

template<typename T>
std::string getTypeAsString() {
  return abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);
}


// following is workaround for some compilers (GCC 4.9.2)
// equivalent to
// template<class...>
// using void_t = void;

template<class...> struct voider { using type = void; };
template<class... Args> using void_t = typename voider<Args...>::type;


template<class T, class = void>
struct has_element_type {
    static constexpr auto value = false;
};

template<class T>
struct has_element_type<T, void_t<typename T::element_type>> {
    static constexpr auto value = true;
};

template<typename T>
typename std::remove_pointer<T>::type
dereference(T* t)
{
    return *t;
}

template<typename T>
typename std::enable_if<!has_element_type<T>::value, T&>::type
dereference(T& t)
{
    return t;
}

// element_type is provided by smart pointers unique_ptr, shared_ptr
template<typename T>
typename std::enable_if<has_element_type<T>::value, typename T::element_type&>::type
dereference(const T& t)
{
    return *t;
}


class execute_on_destroy {
    std::function<void(void)> fct;
public:
    execute_on_destroy(std::function<void(void)> function) : fct(function) {}
    ~execute_on_destroy() {
        fct();
    }
};

}} // namespace ant::std_ext
