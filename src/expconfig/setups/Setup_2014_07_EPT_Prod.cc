#include "Setup_2014_EPT.h"

#include "base/Logger.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_07_EPT_Prod : public Setup_2014_EPT
{
public:

    Setup_2014_07_EPT_Prod(const std::string& name, SetupOptPtr opt)
        : Setup_2014_EPT(name, opt)
    {
        if(!IsFlagSet("IncludeBadElements")) {
            // empty elements
            IgnoreDetectorChannels(Detector_t::Type_t::CB, {203,265,267,479,549,565,586,607,677});
            // inverted timewalk spectrum
            IgnoreDetectorChannels(Detector_t::Type_t::CB, {17,582,672,678,696});
            // misc reasons (e.g. bad timewalk spectrum)
            IgnoreDetectorChannels(Detector_t::Type_t::CB, {41,125,418,547,602});

            // no pi0 peak visible in 1CB 1TAPS hist
            IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {2,73,127,137,138,145,200,218,283,291,356,357,364});

            // no nice timing peak
            IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {6, 64, 128, 192, 287, 321, 337, 349});
        }
        else {
            LOG(WARNING) << "DO NOT ignore bad elements";
        }
    }


    bool Matches(const THeaderInfo& header) const override {
        if(!Setup_2014_EPT::Matches(header))
            return false;
        if(!std_ext::time_between(header.Timestamp, "2014-07-29", "2014-08-25"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_07_EPT_Prod)

}}} // namespace ant::expconfig::setup
