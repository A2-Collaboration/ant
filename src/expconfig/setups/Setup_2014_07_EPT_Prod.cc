#include "Setup_2014_EPT.h"

#include "base/Logger.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for the July 2014 End Point Tagger beam time
 * @see https://wwwa2.kph.uni-mainz.de/intern/daqwiki/analysis/beamtimes/2014-07-29
 */
class Setup_2014_07_EPT_Prod : public Setup_2014_EPT
{
public:

    Setup_2014_07_EPT_Prod(const std::string& name, OptionsPtr opt)
        : Setup_2014_EPT(name, opt)
    {
        SetTimeRange("2014-07-29", "2014-08-25");

        // CB
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {203,265,267,479,549,565,607,677});
        CB->SetElementFlag(Detector_t::ElementFlag_t::BadTDC, {623,662,57,59,162,582,586,672,696});
        CB->SetElementFlag(Detector_t::ElementFlag_t::NoCalibFill, {678,17,557});
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken,  {554}); // not present over full beamtime

        // no Pi0 peak even in MC
        // interestingly, those two elements are opposite to each other
        // (maybe some weirdness in detector model/additional material?)
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::NoCalibFill, {128,347});

        // no Pi0 peak on Data
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::NoCalibFill, {219,375});
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::NoCalibUseDefault, {14});


        // no nice timing peak or very low number of entries
        TAPSVeto->SetElementFlag(Detector_t::ElementFlag_t::Broken, {6, 192, 287, 321, 337, 349});
    }


    triggersimu_config_t GetTriggerSimuConfig() const override
    {
        auto conf = Setup_2014_EPT::GetTriggerSimuConfig();
        // from https://github.com/padlarson/a2GoAT/blob/AdlarsonAnalysis/src/AdlarsonPhysics.cc#L4139
        std_ext::insertRange(conf.CBESum_MissingElements, 352, 415);
        return conf;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_07_EPT_Prod)

}}} // namespace ant::expconfig::setup
