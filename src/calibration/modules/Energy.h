#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"

#include "tree/TID.h" // for TKeyValue, TID

#include <memory>

class TH1;

namespace ant {

namespace calibration {

class DataManager;

namespace gui {
class PeakingFitFunction;
}


class Energy :
        public Calibration::Module, // this makes this module abstract
        public ReconstructHook::DetectorReadHits
{

public:
    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits) override;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;

protected:
    Energy(Detector_t::Type_t detectorType,
           std::shared_ptr<DataManager> calmgr,
           Calibration::Converter::ptr_t converter,
           double defaultPedestal,
           double defaultGain,
           double defaultThreshold,
           double defaultRelativeGain,
           Channel_t::Type_t channelType = Channel_t::Type_t::Integral
           );
    virtual ~Energy();

    /**
     * @brief The CalibType struct stores the data
     * for the Updateable interface and for the GUI
     */
    struct CalibType
    {
        // see also implementation of Get method
        const double        DefaultValue;
        std::vector<double> DefaultValues; // if empty, channel-independent DefaultValue is used
        std::vector<double> Values;        // if empty, channel-dependent DefaultValues[ch] is used
        const std::string   Name;
        const std::string   HistogramName;

        double Get(unsigned channel) const;

        CalibType(double defaultValue, const std::string& name) :
            CalibType(defaultValue, name, name)
        {}
        CalibType(double defaultValue, const std::string& name, const std::string& histname) :
            DefaultValue(defaultValue),
            Values(),
            Name(name),
            HistogramName(histname)
        {}
    }; // CalibType

    /**
     * @brief The GUI_CalibType struct is an abstract base class for
     * handling the data storage for a provided CalibType
     */
    struct GUI_CalibType : gui::CalibModule_traits {
        GUI_CalibType(const std::string& basename,
                      CalibType& type,
                      const std::shared_ptr<DataManager>& calmgr,
                      const std::shared_ptr<Detector_t>& detector_,
                      Calibration::AddMode_t mode = Calibration::AddMode_t::StrictRange
                      );

        virtual std::string GetName() const override;
        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual unsigned GetNumberOfChannels() const override;

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;

        virtual void StartSlice(const interval<TID>& range) override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;

    protected:
        CalibType& calibType;
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<Detector_t> detector;

        std::map< unsigned, std::vector<double> > fitParameters;
        std::vector<double> previousValues;

        bool IgnorePreviousFitParameters = false;
        bool UsePreviousSliceParams = false;

        Calibration::AddMode_t addMode;
    }; // GUI_CalibType


    struct GUI_Pedestals : GUI_CalibType {
        GUI_Pedestals(const std::string& basename,
                      CalibType& type,
                      const std::shared_ptr<DataManager>& calmgr,
                      const std::shared_ptr<Detector_t>& detector,
                      std::shared_ptr<gui::PeakingFitFunction> fitfunction);

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
    protected:
        std::shared_ptr<gui::PeakingFitFunction> func;
        gui::CalCanvas* canvas;
        TH1*  h_projection = nullptr;

    }; // GUI_Pedestals

    struct GUI_Banana : GUI_CalibType {
        GUI_Banana(const std::string& basename,
               CalibType& type,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<Detector_t>& detector);

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual void InitGUI(gui::ManagerWindow_traits* window) override;

        virtual DoFitReturn_t DoFit(TH1* hist, unsigned ch) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;

    protected:
        std::shared_ptr<gui::PeakingFitFunction> func;

        gui::CalCanvas* c_fit;
        gui::CalCanvas* c_extra;
        TH1D* h_projection;
        TH2D* banana;

        TH1D* h_relative = nullptr;

        double AutoStopOnChi2 = 6;
    }; // GUI_Banana


    const Detector_t::Type_t DetectorType;
    const Channel_t::Type_t ChannelType; // can be Integral or IntegralShort

    std::shared_ptr<DataManager> calibrationManager;

    const Calibration::Converter::ptr_t Converter;

    CalibType Pedestals;
    CalibType Gains;
    CalibType Thresholds;
    CalibType RelativeGains;

    std::vector<CalibType*> AllCalibrations = {
        std::addressof(Pedestals),
        std::addressof(Gains),
        std::addressof(Thresholds),
        std::addressof(RelativeGains)
    };

};

}}  // namespace ant::calibration
