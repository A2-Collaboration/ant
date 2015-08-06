#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
class PID;
}}

namespace calibration {

class PID_Energy : public Energy
{


public:
    PID_Energy(
            std::shared_ptr<expconfig::detector::PID> pid,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 100,
            double defaultGain = 0.014,
            double defaultThreshold = 0.001,
            double defaultRelativeGain = 1.0);

    virtual ~PID_Energy();


    class ThePhysics : public Physics {

    protected:
        TH2* pedestals = nullptr;

    public:
         ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const Event& event);
        virtual void Finish();
        virtual void ShowResult();
    };


    virtual std::unique_ptr<Physics> GetPhysicsModule();
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override {}

protected:
    std::shared_ptr<expconfig::detector::PID> pid_detector;

};

}} // namespace ant::calibration
