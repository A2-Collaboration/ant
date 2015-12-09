#include "catch.hpp"

#include "analysis/physics/Physics.h"
#include "expconfig/ExpConfig.h"

#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"

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
    WrapTFileOutput outputfile(tmpfile.filename,
                               WrapTFileOutput::mode_t::recreate,
                               true);

    auto physicsnames = calibration->GetPhysicsModules();
    REQUIRE_FALSE(physicsnames.empty());

    for(auto physicsname : calibration->GetPhysicsModules()) {
        REQUIRE_NOTHROW(ant::analysis::PhysicsRegistry::Create(physicsname));
    }


}

void dotest() {

    for(auto setupname : ExpConfig::Setup::GetNames()) {
        auto setup = ExpConfig::Setup::Get(setupname);
        for(auto calibration : setup->GetCalibrations()) {
            INFO("Setup="+setupname+" Calibration="+calibration->GetName());
            checkcalibration(calibration);
        }
    }

}
