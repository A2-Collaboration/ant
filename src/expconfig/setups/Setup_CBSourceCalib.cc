#include "Setup.h"

#include "detectors/CB.h"
#include "detectors/Trigger.h"


using namespace std;

namespace ant {
namespace expconfig {
namespace setup {


class Setup_CBSourceCalib : public Setup
{
public:

    Setup_CBSourceCalib(const std::string& name, OptionsPtr opt) : Setup(name, opt)
    {
        auto cb = make_shared<detector::CB>();
        AddDetector(cb);
    }
    bool Matches(const TID&) const override {
        return false;
    }
    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override {
        Setup::BuildMappings(hit_mappings, scaler_mappings);
        auto& cb_reftiming = detector::Trigger::Reference_CATCH_CBCrate;
        hit_mappings.emplace_back(
                    cb_reftiming.LogicalChannel,
                    cb_reftiming.AcquRawChannel
                    );
    }
};

AUTO_REGISTER_SETUP(Setup_CBSourceCalib)
}
}
}
