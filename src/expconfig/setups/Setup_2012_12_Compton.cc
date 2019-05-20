#include "Setup.h"

#include "detectors/CB.h"
#include "detectors/PID.h"

#include "base/std_ext/math.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2012_12_Compton : public Setup
{
public:

    Setup_2012_12_Compton(const std::string& name, OptionsPtr opt) : Setup(name, opt)
    {
        /// \todo refine time range for this setup describing the 2012-12 Compton beamtime?
        SetTimeRange("2012-12-01", "2012-12-31");

        auto cb = make_shared<detector::CB>();
        AddDetector(cb);

        auto pid = make_shared<detector::PID_2009_07>();
        AddDetector(pid);

        /// \todo add more detectors/calibrations here, @see Setup_2014_EPT.cc
    }

    virtual double GetElectronBeamEnergy() const override {
        /// \todo put correct beam energy here
        return std::numeric_limits<double>::quiet_NaN();
    }


    virtual candidatebuilder_config_t GetCandidateBuilderConfig() const override {
        candidatebuilder_config_t conf;
        conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
        conf.CB_ClusterThreshold = 15;
        conf.TAPS_ClusterThreshold = 20;
        return conf;
    }
};

AUTO_REGISTER_SETUP(Setup_2012_12_Compton)


}}} // namespace ant::expconfig::setup
