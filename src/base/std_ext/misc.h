#pragma once

#include <cxxabi.h>

namespace ant {
namespace std_ext {

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
