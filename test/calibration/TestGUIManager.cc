#include "catch.hpp"

#include "calibration/gui/Manager.h"
#include "calibration/gui/ManagerWindow_traits.h"
#include "calibration/gui/AvgBuffer.h"
#include "calibration/gui/CalCanvas.h"
#include "calibration/Calibration.h"


#include "analysis/physics/Physics.h"
#include "analysis/plot/HistogramFactory.h"

#include "expconfig/ExpConfig.h"
#include "expconfig_helpers.h"

#include "tree/TAntHeader.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"
#include "base/OptionsList.h"

#include "TROOT.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

void dotest();

TEST_CASE("TestCalibrationModules","[calibration]")
{
    test::EnsureSetup();
    dotest();
}

struct ManagerWindowTest : gui::ManagerWindowGUI_traits {

    ManagerWindowTest() {
        // prevent canvases from doint anything...
        gROOT->SetBatch(true);
    }

    virtual ~ManagerWindowTest() {
        for(auto c : calCanvases)
            delete c;
    }

    vector<gui::CalCanvas*> calCanvases;
    virtual gui::CalCanvas* AddCalCanvas(const string& name) override
    {
        /// \todo make those canvases shared_ptr?
        auto c = new gui::CalCanvas(std_ext::formatter() << name << calCanvases.size());
        calCanvases.push_back(c);
        return c;
    }

    vector<string> checkBoxLabels;
    virtual void AddCheckBox(const string& label, bool&) override
    {
        checkBoxLabels.push_back(label);
    }

    vector<string> numberEntryLabels;
    virtual void AddNumberEntry(const string& label, double&) override
    {
        numberEntryLabels.push_back(label);
    }
    virtual void AddNumberEntry(const std::string& label, double,
                                std::function<void(const TGNumberEntry&)>) override
    {
        numberEntryLabels.push_back(label);
    }

    Mode_t mode;
    virtual Mode_t& GetMode() override
    {
        return mode;
    }

    vector<unsigned> set_max_slices;
    vector<unsigned> set_max_channels;
    virtual void SetProgressMax(unsigned slices, unsigned channels) override
    {
        set_max_channels.push_back(channels);
        set_max_slices.push_back(slices);
    }

    vector<unsigned> set_progress_slice;
    vector<unsigned> set_progress_channel;
    virtual void SetProgress(unsigned slice, unsigned channel) override
    {
        set_progress_slice.push_back(slice);
        set_progress_channel.push_back(channel);
    }

    vector<bool> finishModeSet;
    virtual void SetFinishMode(bool flag) override
    {
        finishModeSet.push_back(flag);
    }
};

void run_calibration(std::shared_ptr< Calibration::PhysicsModule> calibration)
{
    auto& setup = ExpConfig::Setup::Get();
    constexpr auto nSlices = 2;
    // create the requested physics classes
    vector<tmpfile_t> tmpfiles;
    for(int slice=0;slice<nSlices;slice++)
    {
        tmpfiles.emplace_back();

        WrapTFileOutput outputfile(tmpfiles.back().filename, true);

        TAntHeader* header = new TAntHeader();
        gDirectory->Add(header);

        auto physicsnames = calibration->GetPhysicsModules();
        REQUIRE_FALSE(physicsnames.empty());

        for(auto physicsname : calibration->GetPhysicsModules()) {
            REQUIRE_NOTHROW(ant::analysis::PhysicsRegistry::Create(physicsname));
        }

        header->CmdLine = "TestGUIManager";
        header->FirstID = TID(slice, 0, {TID::Flags_t::AdHoc});
        header->LastID = TID(slice, 1, {TID::Flags_t::AdHoc});
        header->SetupName = setup.GetName();
    }
    REQUIRE(tmpfiles.size() == nSlices);

    list<unique_ptr<gui::CalibModule_traits> > guis;
    auto Options = make_shared<OptionsList>();
    calibration->GetGUIs(guis,Options);
    REQUIRE_FALSE(guis.empty());


    for(auto& gui : guis) {

        // Different GUIs may create histograms with identical name
        // create context in TFile to cleanup after each GUI
        tmpfile_t tmpfile;
        WrapTFileOutput outputfile(tmpfile.filename, true);

        INFO("GUI="+gui->GetName());
        vector<string> inputfiles;
        for(auto& tmpfile : tmpfiles) {
            inputfiles.emplace_back(tmpfile.filename);
        }

        gui::Manager manager(inputfiles,
                             // averaging is tested somewhere else
                             std_ext::make_unique<gui::AvgBuffer_Sum<TH1>>(),
                             false // do not confirm header mismatch
                             );
        manager.SetModule(move(gui));
        REQUIRE(manager.DoInit(-1));

        ManagerWindowTest window;
        manager.InitGUI(addressof(window));

        // each GUI should at least create one canvas
        REQUIRE(window.calCanvases.size()>0);


        while(true) {
            auto ret = manager.Run();
            if(ret == gui::Manager::RunReturn_t::Exit)
                break;
        }
    }

}

void dotest() {
    SetErrorHandler([] (
                    int level, Bool_t abort, const char *location,
                    const char *msg) {
        // Fitting is not expected to work, suppress output
        if(string(location) == "Fit")
            return;
        DefaultErrorHandler(level, abort, location, msg);
    });


    auto& setup = ExpConfig::Setup::Get();
    unsigned nCalibrations = 0;
    for(auto calibration : setup.GetCalibrations()) {
        cout << calibration->GetName() << endl;
        INFO("Calibration="+calibration->GetName());
        run_calibration(calibration);
        nCalibrations++;
    }
    REQUIRE(nCalibrations==12);
}
