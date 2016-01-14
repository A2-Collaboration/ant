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
        if(!Options->Get<bool>("IncludeBadElements")) {
            // empty elements
            IgnoreDetectorChannels(Detector_t::Type_t::CB, {203,265,267,479,549,565,586,607,677});
            // inverted timewalk spectrum
            IgnoreDetectorChannels(Detector_t::Type_t::CB, {17,582,672,678,696});
            // misc reasons (e.g. bad timewalk spectrum)
            IgnoreDetectorChannels(Detector_t::Type_t::CB, {41,125,418,547,602});
            // missing partially
            IgnoreDetectorChannels(Detector_t::Type_t::CB, {554});

            // no pi0 peak visible in 1CB 1TAPS hist
            IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {
                                       0,1,2,3,64,72,73,74,75,
                                       127,136,137,138,139,143,144,145,146,147,148,149,
                                       200,210,211,216,217,218,219,220,222,273,
                                       283,284,291,295,346,356,357,358,363,364,367,368,
                                       419,429,436,437
                                   });

            // no nice timing peak or very low number of entries
            IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {6, 192, 287, 321, 337, 349});
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
