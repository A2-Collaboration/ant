#include "Setup.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EPT : public Setup
{
public:

    Setup_2014_EPT(const string& name);

    virtual double GetElectronBeamEnergy() const override;

    virtual cluster_thresholds_t GetClusterThresholds() const override;

    bool Matches(const THeaderInfo& header) const override;

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override;

    virtual ExpConfig::Reconstruct::candidatebuilder_config_t GetCandidateBuilderConfig() const override;
};

}}} // namespace ant::expconfig::setup
