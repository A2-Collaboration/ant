#include "Setup_2010_03_Base.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for the September 2010 beam time
 * @see https://wwwa2.kph.uni-mainz.de/intern/daqwiki/analysis/beamtimes/2010-09-13
 */
class Setup_2010_09_Compton : public Setup_2010_03_Base
{
public:

    Setup_2010_09_Compton(const std::string& name, OptionsPtr opt)
        : Setup_2010_03_Base(name, opt)
    {
        SetTimeRange("2010-09-13", "2010-10-04");

        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {518, 540});

        vector<unsigned> switched_off;
        switched_off.reserve(352-272);
        for(unsigned i=272; i<352; ++i) switched_off.push_back(i);

        Tagger->SetElementFlag(Detector_t::ElementFlag_t::Broken, switched_off);
        Tagger->SetElementFlag(Detector_t::ElementFlag_t::Broken, {27});

        TAPSVeto->SetElementFlag(Detector_t::ElementFlag_t::Broken, {263});
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2010_09_Compton)

}}} // namespace ant::expconfig::setup
