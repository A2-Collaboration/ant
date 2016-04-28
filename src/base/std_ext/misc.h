#pragma once

#include <algorithm>
#include <functional>

#include <cxxabi.h>

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

class execute_on_destroy {
    std::function<void(void)> fct;
public:
    execute_on_destroy(std::function<void(void)> function) : fct(function) {}
    ~execute_on_destroy() {
        fct();
    }
};

}} // namespace ant::std_ext
