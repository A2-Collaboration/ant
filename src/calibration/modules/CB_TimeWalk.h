#pragma once

#include "Calibration.h"

class TGraph;

namespace ant{
namespace calibration{

}
}

namespace ant {


namespace expconfig {
namespace detector {
class CB;
}}

namespace calibration {

class DataManager;

namespace gui {
class FitGaus;
}


class CB_TimeWalk : public Calibration::Module
{
public:
    CB_TimeWalk(
            const std::shared_ptr<expconfig::detector::CB>&  pid,
            const std::shared_ptr<DataManager>& calmgr
            );
    virtual ~CB_TimeWalk();

    class ThePhysics : public analysis::Physics {
    protected:
        TH3D* h_timewalk;

    public:
        ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override ;
        virtual void ShowResult() override;
    }; // ThePhysics

    class TheGUI : public gui::Manager_traits {
    protected:
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<expconfig::detector::CB> cb_detector;
        std::shared_ptr<gui::FitGaus> func;

//        gui::CalCanvas* c_singlechannel;
//        gui::CalCanvas* c_result;

//        TH1*  h_projection = nullptr;
//        TGraph* h_result;

        std::map< unsigned, std::vector<double> > fitParameters;


    public:
        TheGUI(const std::string& basename,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<expconfig::detector::CB>& cb
               );

        virtual std::string GetHistogramName() const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI();
        virtual std::list<gui::CalCanvas*> GetCanvases() const;

        virtual void StartRange(const interval<TID>& range);
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel);
        virtual void DisplayFit();
        virtual void StoreFit(unsigned channel);
        virtual bool FinishRange();
        virtual void StoreFinishRange(const interval<TID>& range);
    }; // TheGUI

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override;

    // Updateable_traits interface
    virtual std::vector<std::list<TID> > GetChangePoints() const override;
    virtual void Update(std::size_t index, const TID& id) override;

protected:
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<DataManager> calibrationManager;



};

}}
