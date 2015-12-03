#include "Setup.h"




#include <limits>

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief The Setup_Test class exists to run the test cases below test/
 */
class Setup_Raw :
        public Setup
{
public:

    Setup_Raw(const std::string& name, SetupOptPtr opt) : Setup(name, opt) {

    }

    virtual double GetElectronBeamEnergy() const override {
        return numeric_limits<double>::quiet_NaN();
    }

    bool Matches(const THeaderInfo&) const override {
        // Setup must be manually selected
        // via command line
        return false;
    }

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override
    {

    }
protected:

};

// if you forget this, your setup is never found....
AUTO_REGISTER_SETUP(Setup_Raw)

}}} // namespace ant::expconfig::setup
