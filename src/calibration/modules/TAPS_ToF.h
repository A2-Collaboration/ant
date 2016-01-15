#pragma once

#include "calibration/Calibration.h"

#include "base/std_ext/math.h"
#include "base/Detector_t.h"
#include "base/interval.h"

#include <memory>

class TH1D;

namespace ant {

namespace expconfig {
namespace detector {
struct TAPS;
}}

namespace calibration {

class DataManager;

namespace gui {
class PeakingFitFunction;
}

class TAPS_ToF :
        public Calibration::Module
{

public:
    class TheGUI : public gui::CalibModule_traits
    {
    protected:
        std::shared_ptr<Detector_t> detector;
        std::shared_ptr<DataManager> calmgr;

        std::vector<double> offsets;
        std::map<unsigned,std::vector<double>> fitParams;

        gui::CalCanvas* theCanvas;
        TH1D*  times;
        TH1D*  timePeaks;

        std::shared_ptr<gui::PeakingFitFunction> fitFunction;
        std::vector<double> previousOffsets;

        bool IgnorePreviousFitParameters = false;
    public:
        TheGUI(const std::string& name,
               const std::shared_ptr<Detector_t>& detector_,
               const std::shared_ptr<DataManager>& cDataManager);

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI(gui::ManagerWindow_traits* window) override;

        virtual void StartSlice(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;
    }; // TheGUI

    TAPS_ToF(const std::shared_ptr<expconfig::detector::TAPS>& detector,
             const std::shared_ptr<DataManager>& CalibrationManager);

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis) override;
    virtual std::vector<std::string> GetPhysicsModules() const;

protected:
    std::shared_ptr<expconfig::detector::TAPS> Detector;

    std::shared_ptr<DataManager> calibrationManager;
};

}}  // namespace ant::calibration
