#pragma once

#include "calibration/Calibration.h"

namespace ant {
namespace calibration {

class DataManager;
using detector_ptr_t = std::shared_ptr<const Detector_t>;

/**
 * @brief The CalibType struct stores the data
 * for the Updateable interface and for the GUI
 */
struct CalibType
{
    const std::string   Name;
    const std::string   HistogramName;
    std::vector<double> Values;        // if empty, channel-dependent DefaultValues[ch] is used

    std::function<void(CalibType&)> NotifyLoad; // called if Values were loaded, see Energy::GetLoaders()

    double Get(unsigned channel) const;

    CalibType(const detector_ptr_t& det,
              const std::string& name,
              const std::vector<double>& defaultValues,
              const std::string& histname = "");

    // prevent copy/move
    CalibType(const CalibType&) = delete;
    CalibType& operator=(const CalibType&) = delete;
    CalibType(CalibType&&) = delete;
    CalibType& operator=(CalibType&&) = delete;
private:
    // if size==1, channel-independent DefaultValue is used
    // see also implementation of Get method
    const std::vector<double> DefaultValues;
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
                  const detector_ptr_t& detector_,
                  Calibration::AddMode_t mode = Calibration::AddMode_t::StrictRange
                  );

    virtual std::string GetName() const override;
    virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
    virtual unsigned GetNumberOfChannels() const override;

    virtual void InitGUI(gui::ManagerWindow_traits& window) override;

    virtual void StartSlice(const interval<TID>& range) override;
    virtual void StoreFinishSlice(const interval<TID>& range) override;

protected:
    OptionsPtr options;
    CalibType& calibType;
    const std::shared_ptr<DataManager> calibrationManager;
    const detector_ptr_t detector;

    std::map< unsigned, std::vector<double> > fitParameters;
    std::vector<double> previousValues;

    bool IgnorePreviousFitParameters = false;
    bool UsePreviousSliceParams = false;

    Calibration::AddMode_t addMode;
}; // GUI_CalibType


}} // namespace ant::calibration