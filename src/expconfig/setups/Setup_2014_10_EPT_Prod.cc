#include "Setup_2014_EPT.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_10_EPT_Prod : public Setup_2014_EPT
{
public:

    Setup_2014_10_EPT_Prod(const std::string& name, SetupOptPtr opt)
        : Setup_2014_EPT(name, opt)
    {
        // see https://wwwa2.kph.uni-mainz.de/intern/daqwiki/analysis/beamtimes/2014-10-14
        IgnoreDetectorChannels(Detector_t::Type_t::CB, {17,125,189,265,267,418,456,547,549,557,565,582,586,597,602,672,677,678,696});
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {127,137,138,145,218,346,356,357,364});
        IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {6,287});
    }


    bool Matches(const THeaderInfo& header) const override {
        if(!Setup_2014_EPT::Matches(header))
            return false;
        if(!std_ext::time_between(header.Timestamp, "2014-10-14", "2014-11-03"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_10_EPT_Prod)

}}} // namespace ant::expconfig::setup
