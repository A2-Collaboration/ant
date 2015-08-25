#pragma once
#include <vector>

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
class SlowcontrolVariable_trais {
public:
    virtual T operator() () const =0;
};

template <typename T>
class ReqestableVariableBase: public SlowcontrolRequestable, public SlowcontrolVariable_trais<T> {
};

template <typename T>
class ReqestableVariable: public ReqestableVariableBase<T> {
protected:
    T data;
public:
    T operator ()() const { return data; }
};

template <typename T>
class ReqestableArray: public ReqestableVariableBase<std::vector<T>> {
protected:
    std::vector<T> data;
public:
    std::vector<T> operator ()() const { return data; }
};


class Slowcontrol {
public:
    ReqestableVariable<double> TotalLivetime;
    ReqestableArray<double> TaggerScalers;
    ReqestableArray<double> EPTScalers;
};
}
}
}
