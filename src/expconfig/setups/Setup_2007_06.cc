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
        SetTimeRange("2007-06-06", "2007-06-24");

        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken,     {518, 540});
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken,     {125}); // uncalibrateable

        // Tagger sections were switched off
        vector<unsigned> switched_off;
        switched_off.reserve(352-224);
        for(unsigned i=224; i<352; ++i) switched_off.push_back(i);
        Tagger->SetElementFlag(Detector_t::ElementFlag_t::Broken, switched_off);

        // Noisy channel
        Tagger->SetElementFlag(Detector_t::ElementFlag_t::Broken, {27});

        //TAPS: No Peak in calibration
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::Broken, {55, 62, 63, 121, 127 ,190, 191, 247, 255, 301, 302, 311, 312, 313, 318, 319, 365});

        //TAPS: No Entries?
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::Broken, {109, 119, 120, 173, 383});

        // TAPSVeto: no entries in calibrationn
        TAPSVeto->SetElementFlag(Detector_t::ElementFlag_t::Broken, {27, 28, 55, 56, 62, 63, 109, 119, 120, 121, 127, 162, 173, 190, 191, 247, 253, 255, 301, 302, 311, 312, 313, 318, 319, 365, 375, 383});

    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2007_06)

}}} // namespace ant::expconfig::setup
