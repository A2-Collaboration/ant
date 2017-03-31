#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"
#include "expconfig/setups/Setup.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/Logger.h"
#include <iostream>

using namespace std;
using namespace ant;

/**
 * @brief Empty fake setup that always matches. Used in case no real setup exists.
 */
struct FakeSetup : ant::expconfig::Setup {
    FakeSetup(): Setup("FakeSetup", make_shared<OptionsList>()) {}
    bool Matches(const TID&) const override { return true; }
    std::list<std::shared_ptr<Detector_t> > GetDetectors() const override { return {}; }
};

int main(int argc, char** argv)
{
    SetupLogger();

    if(argc == 2) {

        // Try to find a setup for the file. Create an empty one if nothing was found.
        const auto getUnpacker = [] (const string& filename) {
            try {
                auto unpacker = Unpacker::Get(filename);
                LOG(INFO) << "Found setup: " <<  ExpConfig::Setup::Get().GetName();
                return unpacker;
            } catch(const ant::ExpConfig::Exception&) {

                LOG(INFO) << "No setup found. Using a fake one. Expect errors!";
            }

            auto setup = make_shared<FakeSetup>();
            expconfig::SetupRegistry::AddSetup(setup->GetName(), setup);
            ExpConfig::Setup::SetByName(setup->GetName());

            return Unpacker::Get(filename);

        };

        auto unpacker = getUnpacker(argv[1]);

        auto firstevent = unpacker->NextEvent();
        if(!firstevent) {
            LOG(ERROR) << "Error: Did not find header info in file " << argv[1];
            return EXIT_FAILURE;
        }
        cout << "TID=" << firstevent.Reconstructed().ID << endl;
        for(auto& message : firstevent.Reconstructed().UnpackerMessages) {
            cout << message << endl;
        }
        return EXIT_SUCCESS;
    }

    cout << "Try to figure out what file that is..." << endl;
    cout << "Usage: Ant-info <filename>" << endl;

    return EXIT_FAILURE;
}
