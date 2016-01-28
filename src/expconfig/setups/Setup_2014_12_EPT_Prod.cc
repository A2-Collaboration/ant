#include "Setup_2014_EPT.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_12_EPT_Prod : public Setup_2014_EPT
{
public:

    Setup_2014_12_EPT_Prod(const std::string& name, SetupOptPtr opt)
        : Setup_2014_EPT(name, opt)
    {
        IgnoreDetectorChannel(Detector_t::Type_t::CB,  17); /// odd time [walk]
        IgnoreDetectorChannel(Detector_t::Type_t::CB,  26); /// empty
        IgnoreDetectorChannel(Detector_t::Type_t::CB,  41); /// odd time walk and energy
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 125); /// odd time [walk] and energy
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 265); /// few entries
        IgnoreDetectorChannel(Detector_t::Type_t::CB, 418); /// odd time spectrum
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

        //IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {36,41,195,203,242,243,254,256,288,292,307,337,349,356,                /// few stat
                                                             //128,129,130,287,320,                                                   /// noise
                                                             //192,263,321});                                                         /// empty
    }


    bool Matches(const THeaderInfo& header) const override {
        if(!Setup_2014_EPT::Matches(header))
            return false;
        if(!std_ext::time_between(header.Timestamp, "2014-12-01", "2014-12-22"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_12_EPT_Prod)

}}} // namespace ant::expconfig::setup
