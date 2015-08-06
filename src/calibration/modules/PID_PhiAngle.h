#pragma once

#include "Calibration.h"

namespace ant{
namespace calibration{
class DataManager;
}
}

namespace ant {


namespace expconfig {
namespace detector {
class PID;
}}

namespace calibration {

class PID_PhiAngle : public Calibration::Module
{
public:
    PID_PhiAngle(std::shared_ptr<expconfig::detector::PID>  pid);
    virtual ~PID_PhiAngle();

    class ThePhysics : public Physics {
    protected:
        TH2* pid_cb_phi_corr;
    public:
        ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const Event& event) override;
        virtual void Finish() override ;
        virtual void ShowResult() override;
    };

    virtual std::unique_ptr<Physics> GetPhysicsModule();
    virtual std::list<std::unique_ptr<calibration::gui::Manager_traits>> GetGUIs() override {
        return {};
    }

    // Updateable_traits interface
    virtual std::vector<std::list<TID> > GetChangePoints() const override;
    virtual void Update(std::size_t index, const TID& id) override;



protected:
    std::shared_ptr<DataManager> calibrationManager;
    std::shared_ptr<expconfig::detector::PID> pid_detector;



};

}}
