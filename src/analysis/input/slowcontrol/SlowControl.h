#pragma once
#include <vector>
#include <cstdint>

namespace ant {
namespace analysis {
namespace input {

class SlowControlRequestable {
protected:
    bool requested = false;
public:
    bool isRequested() const;
    void Request();
};

template <typename T>
class ReqestableVariable: public SlowControlRequestable {
protected:
    T data;
public:
    ReqestableVariable(): data({}) {}

    const T& operator()() const { return data; }
          T& operator()() { return data; }
};

template <> ReqestableVariable<double>::ReqestableVariable();
template <> ReqestableVariable<std::int64_t>::ReqestableVariable();


class SlowControl {
public:
    SlowControl() = default;

    ReqestableVariable<double> TotalLivetime;
    ReqestableVariable<double> FaradayCup;

    ReqestableVariable<std::vector<double>> TaggerScalers;
};
}
}
}
