#pragma once
#include <vector>
#include <cstdint>

namespace ant {
namespace analysis {
namespace data {

class SlowcontrolRequestable {
protected:
    bool requested = false;
public:
    bool isRequested() const;
    void Request();
};

template <typename T>
class ReqestableVariable: public SlowcontrolRequestable {
protected:
    T data;
public:
    ReqestableVariable(): data({}) {}

    const T& operator()() const { return data; }
          T& operator()() { return data; }
};

template <> ReqestableVariable<double>::ReqestableVariable();
template <> ReqestableVariable<std::int64_t>::ReqestableVariable();


class Slowcontrol {
public:
    Slowcontrol() = default;

    ReqestableVariable<double> TotalLivetime;
    ReqestableVariable<double> FaradayCup;

    ReqestableVariable<std::vector<double>> TaggerScalers;
};
}
}
}
