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
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {265,549,565,597,677});
        CB->SetElementFlag(Detector_t::ElementFlag_t::BadTDC, {547,662,678,17,59,162,557,582,586,672,696});
        CB->SetElementFlag(Detector_t::ElementFlag_t::NoCalibFill, {17,678});
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::BadTDC, {137});
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::NoCalibFill, {74,76,148,303,347,353});
        TAPSVeto->SetElementFlag(Detector_t::ElementFlag_t::Broken, {6,64,128,192,256,263,287,321,337,349});
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
