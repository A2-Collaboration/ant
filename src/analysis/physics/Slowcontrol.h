#pragma once

#include "tree/TSlowControl.h"

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
    std::list<Variable*> GetRequestedSlowcontrols() const { return requested_slowcontrols; }
};


class Variable {
private:
    const ant::TSlowControl::Key key;

public:
    Variable(Receiver* receiver, const ant::TSlowControl::Key& Key);

    const ant::TSlowControl::Key GetKey() const { return key; }

};


class Livetime: public Variable {
public:
    Livetime(Receiver* receiver);
};
}
}
}
