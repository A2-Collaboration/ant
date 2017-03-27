#pragma once

#include "calibration/Calibration.h"
#include "fitfunctions/FitGaus.h"

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
class CBPeakFunction: public FitGaus
{
public:
    using FitGaus::FitGaus;
    virtual  void SetDefaults(TH1 *hist) override;
};

}

class Time :
        public Calibration::Module,
        public ReconstructHook::DetectorReadHits
{

public:
    class TheGUI : public gui::CalibModule_traits
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

        bool IgnorePreviousFitParameters = false;
        bool SkipEmptyChannels = false;

        double AutoStopOnChi2 = 6.0;
        double AutoStopOnPeakPos = 6;
        double AutoStopAtChannel = -1;
        double HardTimeCut = -1;
        double Rebin = 1;

        bool channelWasEmpty = false;


    public:
        TheGUI(const std::string& name,
               const std::shared_ptr<Detector_t>& theDetector,
               const std::shared_ptr<DataManager>& cDataManager,
               const std::vector<double>& DefaultOffsets,
               const std::vector<double>& Offsets,
               const std::shared_ptr<gui::PeakingFitFunction> fitFunction);

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI(gui::ManagerWindow_traits& window) override;

        virtual void StartSlice(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(const TH1& hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;
    }; // TheGUI

    Time(const std::shared_ptr<Detector_t>& detector,
         const std::shared_ptr<DataManager>& CalibrationManager,
         Calibration::Converter::ptr_t converter,
         double defaultOffset,
         std::shared_ptr<gui::PeakingFitFunction> FitFunction,
         const interval<double>& timeWindow = {-std_ext::inf, std_ext::inf},
         double defaultGain = 1.0 // default gain is 1.0
         );

    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits) override;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;
    void UpdatedTIDFlags(const TID& id) override;

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;

protected:

    std::shared_ptr<Detector_t> Detector;

    std::shared_ptr<DataManager> calibrationManager;

    std::vector<Calibration::Converter::ptr_t> Converters;

    std::vector<interval<double>> TimeWindows;

    std::shared_ptr<gui::PeakingFitFunction> fitFunction;

    std::vector<double> DefaultOffsets;
    std::vector<double> Offsets;

    std::vector<double> DefaultGains;
    std::vector<double> Gains;

    bool IsMC = false;
};

}}  // namespace ant::calibration
