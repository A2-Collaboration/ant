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
        IgnoreDetectorChannels(Detector_t::Type_t::CB,{  17,     /// odd time [walk]
                                                         26,     /// empty
                                                         41,     /// odd time walk and energy
                                                        125,     /// odd time [walk] and energy
                                                        265,     /// few entries
                                                        418,     /// no peak in Energy calibration
                                                        547,     /// few entries
                                                        549,     /// no entries
                                                        557,     /// double peak in time
                                                        565,     /// no entries
                                                        582,     /// odd time [walk]
                                                        586,     /// few entries
                                                        597,     /// no entries
                                                        602,     /// odd time walk
                                                        672,     /// odd time walk
                                                        677,     /// few entries
                                                        678,     /// inverted/odd spectra
                                                        679,     /// no entries
                                                        696  }); /// odd time walk

        IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {36,41,195,203,242,243,254,256,288,292,307,337,349,356,                /// few stat
                                                             128,129,130,287,320,                                                   /// noise
                                                             192,263,321});                                                         /// empty

        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {
                                   128, 347, // broken in MC
                                   1, 74, 150, 365, 369        // no peak
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
