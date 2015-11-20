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


    class ThePhysics : public analysis::Physics {

    protected:
        TH2D* h_pedestals = nullptr;

        struct PerChannel_t {
            TH2D* PedestalTiming = nullptr;
            TH1D* PedestalNoTiming = nullptr;
            TH2D* Banana = nullptr;
            TH2D* BananaRaw = nullptr;
            TH1D* TDCMultiplicity;
            TH1D* QDCMultiplicity;
            //TH3D* BananaTiming = nullptr;
            PerChannel_t(analysis::SmartHistFactory HistFac);
        };

        std::vector<PerChannel_t> h_perChannel;

    public:
         ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;
    }; // ThePhysics


    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis) override;

protected:
    std::shared_ptr<expconfig::detector::PID> pid_detector;

};

}} // namespace ant::calibration
