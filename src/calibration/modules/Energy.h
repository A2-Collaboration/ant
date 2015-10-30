#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"

#include "tree/TDataRecord.h" // for TKeyValue, TID

#include <memory>

class TH1;

namespace ant {

namespace calibration {

class DataManager;

namespace gui {
class FitGausPol0;
}


class Energy :
        public Calibration::Module, // this makes this module abstract
        public ReconstructHook::DetectorReadHits
{

public:
    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits, extrahits_t& extrahits) override;

    // Updateable_traits interface
    virtual std::vector<std::list<TID>> GetChangePoints() const override;
    void Update(std::size_t index, const TID& tid) override;

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
     * @brief NeedsPedestals can be overridden by base class to disable pedestal insertion
     * @return true if pedestal values (no gain applied, unsubtracted) are needed
     */
    virtual bool NeedsPedestals() const { return true; }

    /**
     * @brief The CalibType struct stores the data
     * for the Updateable interface and for the GUI
     */
    struct CalibType
    {
        const double        DefaultValue;
        std::vector<double> Values;
        const std::string   Name;
        bool Extendable;
        CalibType(double defaultValue, const std::string& name, bool extendable = false):
            DefaultValue(defaultValue),
            Values(),
            Name(name),
            Extendable(extendable)
        {}
    }; // CalibType

    /**
     * @brief The GUI_CalibType struct is an abstract base class for
     * handling the data storage for a provided CalibType
     */
    struct GUI_CalibType : gui::Manager_traits {
        GUI_CalibType(const std::string& basename,
                      CalibType& type,
                      const std::shared_ptr<DataManager>& calmgr,
                      const std::shared_ptr<Detector_t>& detector_
                      );

        virtual std::string GetName() const override;
        virtual std::string GetHistogramName() const override;
        virtual unsigned GetNumberOfChannels() const override;

        virtual void StartRange(const interval<TID>& range) override;
        virtual void StoreFinishRange(const interval<TID>& range) override;

        static std::string ConstructName(const std::string& basename, const std::string& type_name) {
            return basename+"/"+type_name;
        }

    protected:
        CalibType& calibType;
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<Detector_t> detector;

        std::map< unsigned, std::vector<double> > fitParameters;
        std::vector<double> previousValues;

    }; // GUI_CalibType


    struct GUI_Pedestals : GUI_CalibType {
        GUI_Pedestals(const std::string& basename,
               CalibType& type,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<Detector_t>& detector);

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel,
                                    const Manager_traits::DoFitOptions_t& options) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishRange() override;
    protected:
        std::shared_ptr<gui::FitGausPol0> func;
        gui::CalCanvas* canvas;
        TH1*  h_projection = nullptr;

    }; // GUI_Pedestals


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
