#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"
#include "base/OptionsList.h"

#include "tree/TID.h" // for TKeyValue, TID

#include <memory>

class TH1;

namespace ant {

namespace calibration {

class DataManager;

namespace gui {
class PeakingFitFunction;
}

namespace gui {
class FitLandauExpo;
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
           const std::shared_ptr<DataManager>& calmgr,
           const Calibration::Converter::ptr_t& converter,
           std::vector<double> defaultPedestals,
           std::vector<double> defaultGains,
           std::vector<double> defaultThresholds,
           std::vector<double> defaultRelativeGains,
           Channel_t::Type_t channelType = Channel_t::Type_t::Integral
           );
    virtual ~Energy();

    /**
     * @brief The CalibType struct stores the data
     * for the Updateable interface and for the GUI
     */
    struct CalibType
    {
        const std::string   Name;
        const std::string   HistogramName;

        // see also implementation of Get method
        std::vector<double> DefaultValues; // if empty, channel-independent DefaultValue is used
        std::vector<double> Values;        // if empty, channel-dependent DefaultValues[ch] is used

        std::function<void(CalibType&)> NotifyLoad; // called if Values were loaded, see Energy::GetLoaders()

        double Get(unsigned channel) const;

        CalibType(const std::string& name, const std::vector<double>& defaultValues, const std::string& histname = "") :
            Name(name),
            HistogramName(histname.empty() ? name : histname),
            DefaultValues(defaultValues),
            Values()
        {}

        // prevent copy/move
        CalibType(const CalibType&) = delete;
        CalibType& operator=(const CalibType&) = delete;
        CalibType(CalibType&&) = delete;
        CalibType& operator=(CalibType&&) = delete;
    }; // CalibType

    /**
     * @brief The GUI_CalibType struct is an abstract base class for
     * handling the data storage for a provided CalibType
     */
    struct GUI_CalibType : gui::CalibModule_traits {
        GUI_CalibType(const std::string& basename,
                      OptionsPtr options,
                      CalibType& type,
                      const std::shared_ptr<DataManager>& calmgr,
                      const std::shared_ptr<const Detector_t>& detector_,
                      Calibration::AddMode_t mode = Calibration::AddMode_t::StrictRange
                      );

        virtual std::string GetName() const override;
        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual unsigned GetNumberOfChannels() const override;

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;

        virtual void StartSlice(const interval<TID>& range) override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;

    protected:
        OptionsPtr options;
        CalibType& calibType;
        const std::shared_ptr<DataManager> calibrationManager;
        const std::shared_ptr<const Detector_t> detector;

        std::map< unsigned, std::vector<double> > fitParameters;
        std::vector<double> previousValues;

        bool IgnorePreviousFitParameters = false;
        bool UsePreviousSliceParams = false;

        Calibration::AddMode_t addMode;
    }; // GUI_CalibType


    struct GUI_Pedestals : GUI_CalibType {
        GUI_Pedestals(const std::string& basename,
                      OptionsPtr options,
                      CalibType& type,
                      const std::shared_ptr<DataManager>& calmgr,
                      const std::shared_ptr<const Detector_t>& detector,
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
                   OptionsPtr options,
                   CalibType& type,
                   const std::shared_ptr<DataManager>& calmgr,
                   const std::shared_ptr<const Detector_t>& detector,
                   const interval<double>& projectionrange,
                   const double proton_peak_mc_pos
                   );

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

        interval<double> projection_range;
        const double proton_peak_mc;

        double AutoStopOnChi2 = 6;
        const std::string full_hist_name;
    }; // GUI_Banana


    /**
     * @brief MIP stands for Minimum Ionizing Particle
     * and uses electron and positron energy deposits in the PID
     */
    struct GUI_MIP : GUI_CalibType {
        GUI_MIP(const std::string& basename,
                OptionsPtr options,
                CalibType& type,
                const std::shared_ptr<DataManager>& calmgr,
                const std::shared_ptr<const Detector_t>& detector,
                const double peak_mc_pos
                );

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual void InitGUI(gui::ManagerWindow_traits* window) override;

        virtual DoFitReturn_t DoFit(TH1* hist, unsigned ch) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;

    protected:
        std::shared_ptr<gui::PeakingFitFunction> func;

        gui::CalCanvas* canvas;
        TH1D* h_projection = nullptr;
        TH1D* h_peaks = nullptr;
        TH1D* h_relative = nullptr;

        const double peak_mc;

        double AutoStopOnChi2 = 6;
        const std::string full_hist_name;
    }; // GUI_MIP


    /**
     * @brief HEP stands for High Energy Protons
     * and uses the high energy proton tail to calibrate the PID energy
     */
    struct GUI_HEP : GUI_CalibType {
        GUI_HEP(const std::string& basename,
                OptionsPtr options,
                CalibType& type,
                const std::shared_ptr<DataManager>& calmgr,
                const std::shared_ptr<Detector_t>& detector,
                const double proton_peak_mc_pos
                );

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual void InitGUI(gui::ManagerWindow_traits* window) override;

        virtual DoFitReturn_t DoFit(TH1* hist, unsigned ch) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;

    protected:
        std::shared_ptr<gui::PeakingFitFunction> func;

        gui::CalCanvas* canvas;
        TH1D* h_projection = nullptr;
        TH1D* h_peaks = nullptr;
        TH1D* h_relative = nullptr;

        const double proton_peak_mc;

        double AutoStopOnChi2 = 6;
        const std::string full_hist_name;
    }; // GUI_HEP


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
