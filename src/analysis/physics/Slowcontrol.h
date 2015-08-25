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
    std::list<Variable*> GetRequestedSlowcontrols() const { return requested_slowcontrols; }
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


class Livetime: public Variable {
public:
    class MyGetter: public DataGetter {
    protected:
        double v = 0;
    public:
        MyGetter() {}
        virtual ~MyGetter() {}
        virtual void Process(double d) override { v= d*2;}
        virtual std::list<ant::TSlowControl::Key> GetRequiredKeys() const override { return {TSlowControl::Key(TSlowControl::Type_t::EpicsScaler, "TRIG:TotalLivetime")};}
        double get() const { return v; }
    };
protected:


    std::shared_ptr<MyGetter> getter = nullptr;

public:
    using Variable::Variable;
    double operator() () const { return getter->get(); }

    std::shared_ptr<DataGetter> Getter() override {
        getter = std::make_shared<Livetime::MyGetter>();
        return getter;
    }
};


class EpicsPV: public Variable {
protected:
    const std::string pv;
public:

    class MyGetter: public DataGetter {
    protected:
        const std::string pv;
        double v = 0;
    public:
        MyGetter(const std::string& pvname):pv(pvname) {}
        virtual ~MyGetter() {}
        virtual void Process(double d) override { v=d;}
        virtual std::list<ant::TSlowControl::Key> GetRequiredKeys() const override { return {TSlowControl::Key(TSlowControl::Type_t::EpicsScaler, pv)};}
        double get() const { return v; }
    };
protected:

    std::shared_ptr<MyGetter> getter = nullptr;

public:
    EpicsPV(Receiver* receiver, const std::string& pvname):Variable(receiver), pv(pvname) {}
    double operator() () const { return getter->get(); }

    std::shared_ptr<DataGetter> Getter() override {
        getter = std::make_shared<MyGetter>(pv);
        return getter;
    }

    std::string GetEpicsPVName() const { return pv; }
};

class TotalLivetime: public EpicsPV {
public:
    TotalLivetime(Receiver* receiver): EpicsPV(receiver,"TRIG:TotalLivetime") {}
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
