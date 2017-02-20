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
        IgnoreDetectorChannel(Detector_t::Type_t::CB,  17); /// odd time [walk]
        IgnoreDetectorChannel(Detector_t::Type_t::CB,  26); /// empty
        IgnoreDetectorChannel(Detector_t::Type_t::CB,  41); /// odd time walk and energy
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 125); /// odd time [walk] and energy
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 265); /// few entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 547); /// few entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 549); /// no entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 557); /// double peak in time
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 565); /// no entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 582); /// odd time [walk]
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 586); /// few entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 597); /// no entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 602); /// odd time walk
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 672); /// odd time walk
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 677); /// few entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 678); /// inverted/odd spectra
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 679); /// no entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 696); /// odd time walk

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
