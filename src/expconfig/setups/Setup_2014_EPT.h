#include "Setup.h"

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Common base class for all Setups of the 2014 Ent Point Tagger beam times
 */
class Setup_2014_EPT : public Setup
{
public:

    Setup_2014_EPT(const std::string& name, OptionsPtr opt);

    virtual double GetElectronBeamEnergy() const override;

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override;

    virtual ExpConfig::Setup::candidatebuilder_config_t GetCandidateBuilderConfig() const override;

    virtual UnpackerA2GeantConfig::promptrandom_config_t GetPromptRandomConfig() const override;
};

}}} // namespace ant::expconfig::setup
