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
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {2, 73, 127, 137, 138, 144, 145, 146, 147, 200, 218, 220, 222, 283, 291, 346, 353, 356, 357, 363, 364, 368, 437});

        // Not enough statistics for runbyrun "-a 20"
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {0,64,72,75,139,217,295,358,419});

        IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {6,64,128,192,256,263,287,321,337,349});
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
