#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
class TAPSVeto;
}}

namespace calibration {

class TAPSVeto_Energy : public Energy
{


public:
    class ThePhysics : public analysis::Physics {

    protected:
        TH2D* h_pedestals = nullptr;
        TH3D* h_bananas = nullptr;

    public:
        ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;
    };

    TAPSVeto_Energy(std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto,
                    std::shared_ptr<DataManager> calmgr,
                    Calibration::Converter::ptr_t converter,
                    double defaultPedestal = 100,
                    double defaultGain = 0.010,
                    double defaultThreshold = 0.1,
                    double defaultRelativeGain = 1.0);

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override;

protected:
    std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto_detector;
};

}} // namespace ant::calibration
