#include "Setup.h"

#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"

#include "base/std_ext/math.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for Rare Pion Test Beamtime in January 2015 & MC studies
 */


class Setup_2015_01_Pion : public Setup
{
public:

    Setup_2015_01_Pion(const std::string& name, OptionsPtr opt) : Setup(name, opt)
    {
        auto cb = make_shared<detector::CB>();
        AddDetector(cb);

        auto pid = make_shared<detector::PID_2014>();
        AddDetector(pid);

        const bool cherenkovInstalled = false;
        auto taps = make_shared<detector::TAPS_2013>(cherenkovInstalled, false);
        AddDetector(taps);
        auto tapsVeto = make_shared<detector::TAPSVeto_2014>(cherenkovInstalled);
        AddDetector(tapsVeto);

    }

    bool Matches(const TID& tid) const override {
        if(!Setup::Matches(tid))
            return false;
        if(!std_ext::time_between(tid.Timestamp, "2015-01-27", "2015-02-01"))
            return false;
        return true;
    }

    virtual double GetElectronBeamEnergy() const {
        return 450.0;
    }
};



AUTO_REGISTER_SETUP(Setup_2015_01_Pion)

}
}
}
