#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"

#include "tree/TDataRecord.h" // for TKeyValue, TID
#include "base/interval.h"

#include "calibration/fitfunctions/FitGaus.h"

#include <memory>
#include <limits>

#include "TH1D.h"
#include "TH2D.h"

class TH1;

namespace ant {
namespace calibration {

class DataManager;

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

        const double defaultOffset;
        std::vector<double> offsets;
        std::map<unsigned,std::vector<double>> fitParams;

        gui::CalCanvas* theCanvas;
        TH1*  times;
        TH1*  timePeaks;
//        TH2*  timesVSchannels;

        std::shared_ptr<gui::FitGaus> fitFunction;
        std::vector<double> previousValues;

        // Manager_traits interface
    public:
        TheGUI(const std::string& name,
               const std::shared_ptr<Detector_t>& theDetector,
               const std::shared_ptr<DataManager>& cDataManager,
               double DefaultOffset,
               const std::vector<double>& Offsets
               );

        virtual std::string GetHistogramName() const override { return GetName()+"/Offsets";}
        virtual unsigned GetNumberOfChannels() const override { return detector->GetNChannels();}
        virtual void InitGUI() override;
        virtual std::list<gui::CalCanvas*> GetCanvases() const override { return {theCanvas};}

        virtual void StartRange(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishRange() override;
        virtual void StoreFinishRange(const interval<TID>& range) override;
    };

    class ThePhysics : public analysis::Physics {
    public:
        using Physics::Physics;

        ThePhysics(const std::string& name, const std::string& histName,
                   const std::shared_ptr<Detector_t>& theDetector);
        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;

    protected:
        TH2D* hTime;
        std::shared_ptr<Detector_t> detector;
    };

    Time(const std::shared_ptr<Detector_t>& detector,
         const std::shared_ptr<DataManager>& CalibrationManager,
         Calibration::Converter::ptr_t converter,
         double defaultOffset,
         const interval<double>& timeWindow = {-std_ext::inf, std_ext::inf},
         double defaultGain = 1.0, // default gain is 1.0
         const std::vector< TKeyValue<double> >& gains = {}
                                                         );

    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;

    // Updateable_traits interface
    virtual std::vector<std::list<TID>> GetChangePoints() const override;
    void Update(std::size_t index, const TID&) override;


    // Physics_traits interface
    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override {
        return std_ext::make_unique<ThePhysics>(GetName(), "Offsets", Detector);
    }

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override {
        guis.emplace_back(std_ext::make_unique<TheGUI>(
                              GetName(),
                              Detector,
                              calibrationManager,
                              DefaultOffset,
                              Offsets
                              ));
    }

protected:

    std::shared_ptr<Detector_t> Detector;

    std::shared_ptr<DataManager> calibrationManager;

    const Calibration::Converter::ptr_t Converter;

    const interval<double> TimeWindow;

    const double DefaultOffset;
    std::vector<double> Offsets;

    const double DefaultGain;
    std::vector<double> Gains;

};

}}  // namespace ant::calibration
