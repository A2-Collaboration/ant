#include "Setup.h"

#include "detectors/CB.h"
#include "detectors/PID.h"

#include "base/std_ext/math.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2007_Base : public Setup
{
public:

    Setup_2007_Base(const std::string& name, OptionsPtr opt) : Setup(name, opt)
    {
        auto cb = make_shared<detector::CB>();
        AddDetector(cb);

        /// \todo add more detectors/calibrations here, @see Setup_2014_EPT.cc
    }

    virtual double GetElectronBeamEnergy() const override {
        return 1557.0;
    }

    bool Matches(const TID& tid) const override
    {
        /// \todo This is for 2007-06 for now. for testing.
        if(!std_ext::time_between(tid.Timestamp, "2007-06-06", "2007-06-24"))
            return false;
        return true;
    }


    virtual ExpConfig::Setup::candidatebuilder_config_t GetCandidateBuilderConfig() const override {
        candidatebuilder_config_t conf;
        conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
        conf.CB_ClusterThreshold = 15;
        conf.TAPS_ClusterThreshold = 20;
        return conf;
    }
};

AUTO_REGISTER_SETUP(Setup_2007_Base)


}}} // namespace ant::expconfig::setup
