#include "catch.hpp"

#include "analysis/physics/Physics.h"
#include "calibration/Calibration.h"
#include "expconfig/ExpConfig.h"
#include "expconfig_helpers.h"

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

void checkcalibration(std::shared_ptr< Calibration::PhysicsModule> calibration) {

    // create the requested physics classes
    tmpfile_t tmpfile;
    WrapTFileOutput outputfile(tmpfile.filename,
                               WrapTFileOutput::mode_t::recreate,
                               true);

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
    auto setup = ExpConfig::Setup::GetLastFound();
    unsigned nCalibrations = 0;
    for(auto calibration : setup->GetCalibrations()) {
        INFO("Calibration="+calibration->GetName());
        checkcalibration(calibration);
        nCalibrations++;
    }
    REQUIRE(nCalibrations==12);
}
