#pragma once

#include "Energy.h"

class TH1D;

namespace ant {

class TH2CB;

namespace expconfig {
namespace detector {
struct CB;
}}

namespace calibration {

namespace gui {
class FitGausexpo;
}

class CB_SourceCalib : public Calibration::PhysicsModule,
                       public ReconstructHook::DetectorReadHits

{

public:
    virtual void ApplyTo(const readhits_t &hits) override;


    //GUI
    class TheGUI : public gui::CalibModule_traits {
    protected:
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<expconfig::detector::CB> cb_detector;
        std::shared_ptr<gui::FitGausexpo> func;

        gui::CalCanvas* canvas;
        TH1* h_projection = nullptr;
        TH1D* AmBe_peaks = nullptr;
        TH1D* h_peaks = nullptr;
        TH2CB* AmBe_peaks_cb = nullptr;


        std::map< unsigned, std::vector<double> > fitParameters;

        double AutoStopOnChi2 = 6;
        char Histname[100];

    public:
        TheGUI(const std::string& basename,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<expconfig::detector::CB>& cb);
        virtual ~TheGUI();

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI(gui::ManagerWindow_traits& window) override;

        virtual void StartSlice(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(const TH1& hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;
        virtual std::string GetName() const override;

    };

    CB_SourceCalib(
            std::shared_ptr<expconfig::detector::CB> cb,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter
            );

    virtual ~CB_SourceCalib();


    virtual void GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;



protected:
     std::shared_ptr<expconfig::detector::CB> cb_detector;
     std::shared_ptr<DataManager> calibrationManager;
     const Calibration::Converter::ptr_t Converter;
};

}}
