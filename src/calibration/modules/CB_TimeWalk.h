#pragma once

#include "Calibration.h"

class TGraph;
class TH1D;
class TH2D;
class TF1;

namespace ant {

namespace expconfig {
namespace detector {
struct CB;
}} // namespace expconfig::detector

namespace calibration {

class DataManager;

namespace gui {
class FitTimewalk;
}


class CB_TimeWalk :
        public Calibration::Module,
        public ReconstructHook::ClusterHits
{
public:
    CB_TimeWalk(
            const std::shared_ptr<expconfig::detector::CB>& cb,
            const std::shared_ptr<DataManager>& calmgr,
            const interval<double>& timeWindow,
            double badTDC_EnergyThreshold = std_ext::inf
            );
    virtual ~CB_TimeWalk();

    class TheGUI : public gui::CalibModule_traits {
    protected:
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<expconfig::detector::CB> cb_detector;
        std::vector< std::shared_ptr<gui::FitTimewalk> >& timewalks;
        std::shared_ptr<gui::FitTimewalk> last_timewalk;

        gui::CalCanvas* c_fit;
        gui::CalCanvas* c_extra;
        gui::CalCanvas* c_extra1;
        gui::CalCanvas* c_extra2;
        TH1D* means;
        TH1D* amplitudes;
        TH1D* sigmas;
        TH2D* proj;

        double slicesY_entryCut = 5;
        TF1* slicesY_gaus = nullptr;

    public:
        TheGUI(const std::string& basename,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<expconfig::detector::CB>& cb,
               std::vector<std::shared_ptr<gui::FitTimewalk> >& timewalks_);

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI(gui::ManagerWindow_traits* window) override;


        virtual void StartSlice(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned ch) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;
    }; // TheGUI

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;
    virtual void ApplyTo(clusterhits_t& sorted_clusterhits) override;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;
    void UpdatedTIDFlags(const TID& id) override;


protected:
    std::vector< std::shared_ptr<gui::FitTimewalk> > timewalks;

    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<DataManager> calibrationManager;

    const interval<double> TimeWindow;
    bool IsMC = false;
    const double BadTDC_EnergyThreshold;

};

}}
