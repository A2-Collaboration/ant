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
        // empty elements
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {203,265,267,479,549,565,586,607,677});
        // inverted timewalk spectrum
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {17,582,672,678,696});
        // misc reasons (e.g. bad timewalk spectrum)
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {41,125,547,602});
        // missing partially
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {554});

        // no Pi0 peak even in MC
        // interestingly, those two elements are opposite to each other
        // (maybe some weirdness in detector model/additional material?)
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::Broken, {128,347});

        // no Pi0 peak on Data
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::Broken, {3,219,375});


        // no nice timing peak or very low number of entries
        TAPSVeto->SetElementFlag(Detector_t::ElementFlag_t::Broken, {6, 192, 287, 321, 337, 349});
    }


    bool Matches(const TID& tid) const override {
        if(!std_ext::time_between(tid.Timestamp, "2014-07-29", "2014-08-25"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_07_EPT_Prod)

}}} // namespace ant::expconfig::setup
