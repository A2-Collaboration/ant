#include "catch.hpp"

#include "calibration/gui/Manager.h"
#include "calibration/gui/ManagerWindow_traits.h"
#include "calibration/gui/AvgBuffer.h"
#include "analysis/physics/Physics.h"
#include "calibration/Calibration.h"
#include "expconfig/ExpConfig.h"
#include "expconfig_helpers.h"

#include "tree/TAntHeader.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"
#include "base/OptionsList.h"

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

    virtual gui::CalCanvas* AddCalCanvas(const string& name) override
    {
    }
    virtual void AddCheckBox(const string& label, bool& flag) override
    {
    }
    virtual void AddNumberEntry(const string& label, double& number) override
    {
    }

    virtual Mode_t& GetMode() override
    {
    }
    virtual void SetProgressMax(unsigned slices, unsigned channels) override
    {
    }
    virtual void SetProgress(unsigned slice, unsigned channel) override
    {
    }
    virtual void SetFinishMode(bool flag) override
    {
    }
};

void run_calibration(std::shared_ptr< Calibration::PhysicsModule> calibration)
{
    auto setup = ExpConfig::Setup::GetLastFound();
    constexpr auto nSlices = 10;
    // create the requested physics classes
    vector<tmpfile_t> tmpfiles;
    for(int slice=0;slice<nSlices;slice++)
    {
        tmpfiles.emplace_back();

        WrapTFileOutput outputfile(tmpfiles.back().filename,
                                   WrapTFileOutput::mode_t::recreate,
                                   true);

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
        header->SetupName = setup->GetName();
    }
    REQUIRE(tmpfiles.size() == nSlices);

    list<unique_ptr<gui::CalibModule_traits> > guis;
    auto Options = make_shared<OptionsList>();
    calibration->GetGUIs(guis,Options);
    REQUIRE_FALSE(guis.empty());

    for(auto& gui : guis) {
        INFO("GUI="+gui->GetName());
        vector<string> inputfiles;
        for(auto& tmpfile : tmpfiles) {
            inputfiles.emplace_back(tmpfile.filename);
        }

        gui::Manager manager(inputfiles,
                             std_ext::make_unique<gui::AvgBuffer_SavitzkyGolay>(3, 0),
                             false // do not confirm header mismatch
                             );
        manager.SetModule(move(gui));


    }

}

void dotest() {
    auto setup = ExpConfig::Setup::GetLastFound();
    REQUIRE(setup != nullptr);
    unsigned nCalibrations = 0;
    for(auto calibration : setup->GetCalibrations()) {
        if(calibration->GetName() != "CB_Time")
            continue;
        INFO("Calibration="+calibration->GetName());
        run_calibration(calibration);
        nCalibrations++;
    }
//    REQUIRE(nCalibrations==12);
}
