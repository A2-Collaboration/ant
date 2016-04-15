#include "Setup_2014_EPT.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for the October 2014 End Point Tagger beam time
 * @see https://wwwa2.kph.uni-mainz.de/intern/daqwiki/analysis/beamtimes/2014-10-14
 */
class Setup_2014_10_EPT_Prod : public Setup_2014_EPT
{
public:

    Setup_2014_10_EPT_Prod(const std::string& name, OptionsPtr opt)
        : Setup_2014_EPT(name, opt)
    {
        // see https://wwwa2.kph.uni-mainz.de/intern/daqwiki/analysis/beamtimes/2014-10-14
        IgnoreDetectorChannels(Detector_t::Type_t::CB, {17,125,189,265,267,418,456,547,549,557,565,582,586,597,602,672,677,678,696});
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {2, 73, 127, 137, 138, 144, 145, 146, 147, 200, 218, 220, 222, 283, 291, 346, 353, 356, 357, 363, 364, 368, 437});

        // Not enough statistics for runbyrun "-a 20"
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {0,3,64,72,74,75,139,149,210,217,273,284,295,347,358,362,367,419});

        IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {6,64,128,192,256,263,287,321,337,349});
    }


    bool Matches(const TID& tid) const override {
        if(!std_ext::time_between(tid.Timestamp, "2014-10-14", "2014-11-03"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_10_EPT_Prod)

}}} // namespace ant::expconfig::setup
