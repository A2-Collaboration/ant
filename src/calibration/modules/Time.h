#pragma once

#include "calibration/Calibration.h"

#include "base/std_ext/math.h"
#include "base/Detector_t.h"
#include "base/interval.h"

#include <memory>

class TH1;

namespace ant {
namespace calibration {

class DataManager;

namespace gui {
class PeakingFitFunction;
}

class Time :
        public Calibration::Module,
        public ReconstructHook::DetectorReadHits
{

public:
    class TheGUI : public gui::Manager_traits
    {
    protected:
        std::shared_ptr<Detector_t> detector;
        std::shared_ptr<DataManager> calmgr;

        const std::vector<double> defaultOffsets;
        std::vector<double> offsets;
        std::map<unsigned,std::vector<double>> fitParams;

        gui::CalCanvas* theCanvas;
        TH1*  times;
        TH1*  timePeaks;

        std::shared_ptr<gui::PeakingFitFunction> fitFunction;
        std::vector<double> previousOffsets;

    public:
        TheGUI(const std::string& name,
               const std::shared_ptr<Detector_t>& theDetector,
               const std::shared_ptr<DataManager>& cDataManager,
               const std::vector<double>& DefaultOffsets,
               const std::vector<double>& Offsets,
               const std::shared_ptr<gui::PeakingFitFunction> fitFunction);

        virtual std::string GetHistogramName() const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI(gui::ManagerWindow_traits* window) override;

        virtual void StartRange(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel,
                                    const Manager_traits::DoFitOptions_t& options) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishRange() override;
        virtual void StoreFinishRange(const interval<TID>& range) override;
    }; // TheGUI

    class ThePhysics : public analysis::Physics {
    public:
        ThePhysics(const std::string& name, const std::string& histName,
                   const std::shared_ptr<Detector_t>& theDetector);
        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void ShowResult() override;

    protected:
        TH2D* hTime;
        std::shared_ptr<Detector_t> detector;
        bool isTagger;
    }; // ThePhysics

    Time(const std::shared_ptr<Detector_t>& detector,
         const std::shared_ptr<DataManager>& CalibrationManager,
         Calibration::Converter::ptr_t converter,
         double defaultOffset,
         std::shared_ptr<gui::PeakingFitFunction> FitFunction,
         const interval<double>& timeWindow = {-std_ext::inf, std_ext::inf},
         double defaultGain = 1.0 // default gain is 1.0
         );

    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;

    // Updateable_traits interface
    virtual std::vector<std::list<TID>> GetChangePoints() const override;
    void Update(std::size_t index, const TID&) override;
    void UpdatedTIDFlags(const TID& id) override;


    // Physics_traits interface
    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override;

protected:

    std::shared_ptr<Detector_t> Detector;

    std::shared_ptr<DataManager> calibrationManager;

    const Calibration::Converter::ptr_t Converter;

    std::vector<interval<double>> TimeWindows;

    std::shared_ptr<gui::PeakingFitFunction> fitFunction;

    std::vector<double> DefaultOffsets;
    std::vector<double> Offsets;

    std::vector<double> DefaultGains;
    std::vector<double> Gains;

    bool IsMC = false;
};

}}  // namespace ant::calibration
