#include "catch.hpp"
#include "expconfig_helpers.h"

#include "analysis/physics/Physics.h"
#include "expconfig/ExpConfig.h"

#include "analysis/utils/Uncertainties.h"

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
    test::EnsureSetup();
    dotest();
}

bool histogram_overwrite_detected = false;
bool duplicate_mkdir_detected = false;


void dotest() {

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
    for(auto name : PhysicsRegistry::GetList()) {

        // Exclude known physics classes which implement method calls which will segfault when using ROOT6
        /// \todo Check if the bug reported <a href="https://sft.its.cern.ch/jira/browse/ROOT-9450">here</a> is fixed and remove exclusion
        /// \bug Certain physics classes will probably crash due to <a href="https://sft.its.cern.ch/jira/browse/ROOT-9450">ROOT-9450</a>
        #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
        for (const string& physics : {"EtapDalitz", "EtapOmegaG", "Omega_EpEm", "sigmaPlus", "triplePi0"})
            if (!physics.compare(name))
                continue;
        #endif

        cout << "Running physics class " << name << endl;
        // some errors only appear when some outfiles are present
        tmpfile_t tmpfile;
        auto outfile = std_ext::make_unique<WrapTFileOutput>(tmpfile.filename, true);

        histogram_overwrite_detected = false;
        duplicate_mkdir_detected = false;
        INFO(name);
        try {
            PhysicsRegistry::Create(name);
            REQUIRE_FALSE(histogram_overwrite_detected);
            REQUIRE_FALSE(duplicate_mkdir_detected);
            auto objects = outfile->GetListOf<TNamed>();
            if(objects.size()>1) {
                for(auto o : objects)
                    cout << "Found object " << o->GetName() << endl;
            }
            REQUIRE(objects.size() == 1);
            REQUIRE(objects.front()->GetName() == name);
        }
        catch(PhysicsRegistry::Exception& e) {
            FAIL(string("Physics Registry error: ")+e.what());
        }
        catch(WrapTFile::Exception&) {
            // ignore silently if Physics classes can't load some files...
            continue;
        }
        catch(ExpConfig::ExceptionNoDetector&) {
            // ignore silently if test setup did not provide detector
            continue;
        }
        catch(utils::UncertaintyModel::Exception&) {
            // ignore silently if class cannot load uncertainty model
            continue;
        }
        catch(Physics::ExceptionOptionNeeded&) {
            // ignore silently if class needs user option
            continue;
        }
        catch(const std::exception& e) {
            FAIL(string("Unexpected exception: ")+e.what());
        }
        catch(...) {
            FAIL("Something weird was thrown.");
        }
        // write the file
        REQUIRE_NOTHROW(outfile = nullptr);

    }


}
