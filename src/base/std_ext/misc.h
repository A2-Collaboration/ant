#pragma once

#include <algorithm>
#include <functional>

#include <cxxabi.h>

namespace ant {
namespace std_ext {

using namespace std;


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
