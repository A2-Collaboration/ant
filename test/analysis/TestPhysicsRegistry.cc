#include "catch.hpp"

#include "analysis/physics/Physics.h"
#include "base/OptionsList.h"

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

void dotest() {
    SetErrorHandler([] (
                    int level, Bool_t abort, const char *location,
                    const char *msg) {
        if(string(location) == "TDirectory::Append")
            histogram_overwrite_detected = true;
        DefaultErrorHandler(level, abort, location, msg);
    });
    std::shared_ptr<OptionsList> popts = make_shared<OptionsList>();
    for(auto name : PhysicsRegistry::GetList()) {
        histogram_overwrite_detected = false;
        INFO(name);
        REQUIRE_NOTHROW(PhysicsRegistry::Create(name, popts));
        REQUIRE_FALSE(histogram_overwrite_detected);
    }
}
