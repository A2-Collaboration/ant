#pragma once

#include "tree/TSlowControl.h"

#include <list>
#include <set>
#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {

class Variable;

class Receiver {
    friend class Variable;

private:
    std::list<Variable*> requested_slowcontrols;

protected:
    void RequestSlowcontrol(Variable* var);

public:
    std::list<Variable*> GetRequestedSlowcontrols() const;
};


class DataGetter {
public:
    virtual void Process(double d) =0;
    virtual std::list<ant::TSlowControl::Key> GetRequiredKeys() const =0;
    virtual ~DataGetter() {}
};

class Variable {
public:
    Variable(Receiver* receiver);
    virtual std::shared_ptr<DataGetter> Getter() =0;
};

class Distributor {
public:
    std::list<std::shared_ptr<DataGetter>> getters;
    std::set<ant::TSlowControl::Key> requestedKeys;
    void Register(Receiver& rec);
    void Process(double d);
};
}
}
}
