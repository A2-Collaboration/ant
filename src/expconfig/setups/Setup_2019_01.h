#include "Setup.h"

#include "detectors/Trigger.h"
#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"
#include "detectors/Tagger.h"
#include "detectors/Cherenkov.h"

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief
 */
class Setup_2019_01 : public Setup
{
protected:
    const bool MCTaggerHits;
    const bool cherenkovInstalled;
    const bool pizzaInstalled;
    const std::shared_ptr<detector::Trigger_2014> Trigger;
    const std::shared_ptr<detector::Tagger_2019_01> Tagger;
    const std::shared_ptr<detector::CB> CB;
    const std::shared_ptr<detector::PID_2014> PID;
    const std::shared_ptr<detector::TAPS_2013_11> TAPS;
    const std::shared_ptr<detector::TAPSVeto_2014> TAPSVeto;

public:

    Setup_2019_01(const std::string& name, OptionsPtr opt);

    virtual double GetElectronBeamEnergy() const override;

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override;

    virtual candidatebuilder_config_t GetCandidateBuilderConfig() const override;

    virtual triggersimu_config_t GetTriggerSimuConfig() const override;

    virtual target_properties_t GetTargetProperties() const override;

    virtual UnpackerA2GeantConfig::promptrandom_config_t GetPromptRandomConfig() const override;

    virtual bool Matches(const TID &tid) const override;
};

}}} // namespace ant::expconfig::setup
