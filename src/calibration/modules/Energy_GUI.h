#pragma once

#include "CalibType.h"

class TF1;

namespace ant {
namespace calibration {

class DataManager;

namespace gui {
class PeakingFitFunction;
}

namespace gui {
class FitLandauExpo;
}

namespace gui {
class FitVetoBand;
}

namespace energy {

struct GUI_Pedestals : GUI_CalibType {
    GUI_Pedestals(const std::string& basename,
                  OptionsPtr options,
                  CalibType& type,
                  const std::shared_ptr<DataManager>& calmgr,
                  const detector_ptr_t& detector,
                  std::shared_ptr<gui::PeakingFitFunction> fitfunction);

    virtual void InitGUI(gui::ManagerWindow_traits& window) override;
    virtual DoFitReturn_t DoFit(const TH1& hist, unsigned channel) override;
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
               const detector_ptr_t& detector,
               const interval<double>& projectionrange,
               const double proton_peak_mc_pos
               );

    virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
    virtual void InitGUI(gui::ManagerWindow_traits& window) override;

    virtual DoFitReturn_t DoFit(const TH1& hist, unsigned ch) override;
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
            const detector_ptr_t& detector,
            const double peak_mc_pos
            );

    virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
    virtual void InitGUI(gui::ManagerWindow_traits& window) override;

    virtual DoFitReturn_t DoFit(const TH1& hist, unsigned ch) override;
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
            const detector_ptr_t& detector,
            const double proton_peak_mc_pos
            );

    virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
    virtual void InitGUI(gui::ManagerWindow_traits& window) override;

    virtual DoFitReturn_t DoFit(const TH1& hist, unsigned ch) override;
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

    double AutoStopOnChi2 = 16;
    const std::string full_hist_name;
}; // GUI_HEP


/**
 * @brief fits slices of the PID energy along the energy axis
 * to fit the shape of the PID banana
 */
struct GUI_BananaSlices : GUI_CalibType {
    GUI_BananaSlices(const std::string& basename,
                     OptionsPtr options,
                     CalibType& type,
                     const std::shared_ptr<DataManager>& calmgr,
                     const detector_ptr_t& detector,
                     const interval<double>& fitrange
                     );

    virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
    virtual void InitGUI(gui::ManagerWindow_traits& window) override;

    virtual DoFitReturn_t DoFit(const TH1& hist, unsigned ch) override;
    virtual void DisplayFit() override;
    virtual void StoreFit(unsigned channel) override;
    virtual bool FinishSlice() override;

protected:
    std::shared_ptr<gui::FitVetoBand> func;

    gui::CalCanvas* c_fit;
    gui::CalCanvas* c_extra;
    TH1D* means;
    TH2D* proj;

    TF1* slicesY_gaus = nullptr;
    double AutoStopOnChi2 = 6;
    double slicesY_entryCut = 5;
    double slicesY_IQRFactor_lo = 1;
    double slicesY_IQRFactor_hi = 3;

    //gui::CalCanvas* canvas;
    TH1D* h_projection = nullptr;
    TH1D* h_vals = nullptr;
    TH1D* h_relative = nullptr;

    interval<double> fit_range;

    const std::string full_hist_name;
}; // GUI_BananaSlices

}}} // namespace ant::calibration::energy
