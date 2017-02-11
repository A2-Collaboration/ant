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
        IgnoreDetectorChannels(Detector_t::Type_t::CB, {17,125,189,265,267,456,547,549,557,565,582,586,597,602,672,677,678,696});
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {1,76, 128,149,347,353});
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
