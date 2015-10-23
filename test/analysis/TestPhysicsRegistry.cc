#include "catch.hpp"

#include "analysis/physics/Physics.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"

#include "TError.h"

#include <memory>

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest();

TEST_CASE("PhysicsRegistry: Create all physics classes", "[analysis]") {
    dotest();
}

bool histogram_overwrite_detected = false;
bool duplicate_mkdir_detected = false;


void dotest() {
    // some errors only appear when some outfiles are present
    tmpfile_t tmpfile;
    auto outfile = std_ext::make_unique<WrapTFileOutput>(tmpfile.filename,
                            WrapTFileOutput::mode_t::recreate,
                            true);
    // overwrite ROOT's error handler to detect some warnings
    SetErrorHandler([] (
                    int level, Bool_t abort, const char *location,
                    const char *msg) {
        // those tests are specific enough...
        if(string(location) == "TDirectory::Append")
            histogram_overwrite_detected = true;
        if(string(location) == "TDirectoryFile::Append")
            histogram_overwrite_detected = true;
        if(string(location) == "TDirectoryFile::mkdir")
            duplicate_mkdir_detected = true;
        DefaultErrorHandler(level, abort, location, msg);
    });

    // create all available physics classes
    std::shared_ptr<OptionsList> popts = make_shared<OptionsList>();
    for(auto name : PhysicsRegistry::GetList()) {
        histogram_overwrite_detected = false;
        duplicate_mkdir_detected = false;
        INFO(name);
        REQUIRE_NOTHROW(PhysicsRegistry::Create(name, popts));
        REQUIRE_FALSE(histogram_overwrite_detected);
        REQUIRE_FALSE(duplicate_mkdir_detected);
    }

    // write the file
    REQUIRE_NOTHROW(outfile = nullptr);
}
