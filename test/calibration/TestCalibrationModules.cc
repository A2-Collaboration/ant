#include "catch.hpp"

#include "analysis/physics/Physics.h"
#include "expconfig/ExpConfig.h"

#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"
#include "base/OptionsList.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

void dotest();

TEST_CASE("TestCalibrationModules","[calibration]")
{
    dotest();
}

void checkcalibration(std::shared_ptr< Calibration::PhysicsModule> calibration) {

    // create the requested physics classes
    tmpfile_t tmpfile;
    WrapTFileOutput outputfile(tmpfile.filename, true);

    auto physicsnames = calibration->GetPhysicsModules();
    REQUIRE_FALSE(physicsnames.empty());

    for(auto physicsname : calibration->GetPhysicsModules()) {
        REQUIRE_NOTHROW(ant::analysis::PhysicsRegistry::Create(physicsname));
    }

    list<unique_ptr<gui::CalibModule_traits> > guis;
    auto Options = make_shared<OptionsList>();
    calibration->GetGUIs(guis,Options);
    REQUIRE_FALSE(guis.empty());

    for(const auto& gui : guis) {
        INFO("GUI="+gui->GetName())
        REQUIRE(gui->GetHistogram(outputfile) != nullptr);
    }

}

void dotest() {

    for(auto setupname : ExpConfig::Setup::GetNames()) {
        ExpConfig::Setup::SetByName(setupname);
        auto& setup = ExpConfig::Setup::Get();
        for(auto calibration : setup.GetCalibrations()) {
            INFO("Setup="+setupname+" Calibration="+calibration->GetName());
            checkcalibration(calibration);
        }
    }

}
