#include "Setup_2007_Base.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for the Juni 2007 beam time
 * @see https://wwwa2.kph.uni-mainz.de/intern/daqwiki/analysis/beamtimes/2007-06
 */
class Setup_2007_06 : public Setup_2007_Base
{
public:

    Setup_2007_06(const std::string& name, OptionsPtr opt)
        : Setup_2007_Base(name, opt)
    {
        IgnoreDetectorChannels(Detector_t::Type_t::CB,     {518, 540});
        IgnoreDetectorChannels(Detector_t::Type_t::CB,     {125}); // uncalibrateable

        // Tagger sections were switched off
        vector<unsigned> switched_off;
        switched_off.reserve(352-224);
        for(unsigned i=224; i<352; ++i) switched_off.push_back(i);
        IgnoreDetectorChannels(Detector_t::Type_t::Tagger, switched_off);

        // Noisy channel
        IgnoreDetectorChannel(Detector_t::Type_t::Tagger, 27);

        //TAPS: No Peak in calibration
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {55, 62, 63, 121, 127 ,190, 191, 247, 255, 301, 302, 311, 312, 313, 318, 319, 365});

        //TAPS: No Entries?
        IgnoreDetectorChannels(Detector_t::Type_t::TAPS, {109, 119, 120, 173, 383});

        // TAPSVeto: no entries in calibrationn
        IgnoreDetectorChannels(Detector_t::Type_t::TAPSVeto, {27, 28, 55, 56, 62, 63, 109, 119, 120, 121, 127, 162, 173, 190, 191, 247, 253, 255, 301, 302, 311, 312, 313, 318, 319, 365, 375, 383});

    }


    bool Matches(const TID& tid) const override {
        if(!std_ext::time_between(tid.Timestamp, "2007-06-06", "2007-06-24"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2007_06)

}}} // namespace ant::expconfig::setup
