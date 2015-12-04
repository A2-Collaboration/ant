#include "Setup.h"

#include "base/Logger.h"
#include "base/piecewise_interval.h"

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

    Setup_Raw(const std::string& name, SetupOptPtr opt) : Setup(name, opt)
    {
        ADC_ranges = Options->Get<decltype(ADC_ranges)>("AcquADC");
        if(ADC_ranges.empty())
            LOG(ERROR) << "No Acqu ADC numbers supplied. Won't write anything.";
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
                       std::vector<scaler_mapping_t>&) const override
    {
        // create a 1:1 mapping to LogicalChannels
        // the Detector is set to "Raw"

        for(interval<unsigned> ADC_range : ADC_ranges) {
            if(!ADC_range.IsSane()) {
                LOG(WARNING) << "Skipping invalid Acqu ADC range " << ADC_range;
                continue;
            }
            for(unsigned adc = ADC_range.Start(); adc <= ADC_range.Stop(); adc++) {
                hit_mappings.emplace_back(
                            Detector_t::Type_t::Raw,
                            Channel_t::Type_t::Raw,
                            adc,
                            adc
                            );
            }
            LOG(INFO) << "Added Acqu ADC range " << ADC_range;
        }

    }
protected:

    PiecewiseInterval<unsigned> ADC_ranges;
};

// if you forget this, your setup is never found....
AUTO_REGISTER_SETUP(Setup_Raw)

}}} // namespace ant::expconfig::setup
