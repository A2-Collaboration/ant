#include "Setup_2014_EPT.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for the December 2014 End Point Tagger beam time
 * @see https://wwwa2.kph.uni-mainz.de/intern/daqwiki/analysis/beamtimes/2014-10-14
 */
class Setup_2014_12_EPT_Prod : public Setup_2014_EPT
{
public:

    Setup_2014_12_EPT_Prod(const std::string& name, OptionsPtr opt)
        : Setup_2014_EPT(name, opt)
    {
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {265,549,557,565,597,677});
        CB->SetElementFlag(Detector_t::ElementFlag_t::BadTDC, {662,678,17,59,162,265,418,582,586,672,696});
        CB->SetElementFlag(Detector_t::ElementFlag_t::NoCalib,{678});

        TAPSVeto->SetElementFlag(Detector_t::ElementFlag_t::Broken, {
                                     36,41,195,203,242,243,254,256,288,292,307,337,349,356,  /// few stat
                                     128,129,130,287,320,                                    /// noise
                                     192,263,321
                                 });                                                         /// empty

        TAPS->SetElementFlag(Detector_t::ElementFlag_t::Broken, {
                                   128, 347,                       // broken in MC
                                   1, 74, 150, 365, 369            // no peak
                               });
    }


    bool Matches(const TID& tid) const override {
        if(!std_ext::time_between(tid.Timestamp, "2014-12-01", "2014-12-22"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_12_EPT_Prod)

}}} // namespace ant::expconfig::setup
