#pragma once

#include "analysis/physics/slowcontrol/Slowcontrol.h"

namespace ant {
namespace analysis {
namespace slowcontrol {

/**
 * @brief Generic slowcontrol accessor for EPICS PVs
 */
class EpicsPV: public Variable {
protected:
    const std::string pv;

    class MyGetter: public DataGetter {
    protected:
        const std::string pv;
        double v = 0;
    public:
        MyGetter(const std::string& pvname);
        virtual ~MyGetter() {}
        virtual void Process(double d) override;
        virtual std::list<ant::TSlowControl::Key> GetRequiredKeys() const override;
        double get() const { return v; }
    };

    std::shared_ptr<MyGetter> getter = nullptr;

public:
    EpicsPV(Receiver* receiver, const std::string& pvname);

    double operator() () const { return getter->get(); }

    std::shared_ptr<DataGetter> Getter() override {
        getter = std::make_shared<MyGetter>(pv);
        return getter;
    }

    std::string GetEpicsPVName() const { return pv; }
};

/**
 * @brief Total Livetime accessor via EPICS (TRIG:TotalLivetime)
 */
class TotalLivetime: public EpicsPV {
public:
    TotalLivetime(Receiver* receiver);
};

}
}
}
