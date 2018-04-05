#include "Setup.h"

#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/Tagger.h"
#include "detectors/Trigger.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Common base class for all Setups of the 2014 Ent Point Tagger beam times
 */
class Setup_2007_Base : public Setup
{
protected:
    const bool MCTaggerHits;
    const std::shared_ptr<detector::Trigger_2007> Trigger;
    const std::shared_ptr<detector::Tagger_2007> Tagger;
    const std::shared_ptr<detector::CB> CB;
    const std::shared_ptr<detector::PID_2007> PID;
    const std::shared_ptr<detector::TAPS_2007> TAPS;
    const std::shared_ptr<detector::TAPSVeto_2007> TAPSVeto;
public:

    Setup_2007_Base(const std::string& name, OptionsPtr opt);

    virtual double GetElectronBeamEnergy() const override;

    virtual candidatebuilder_config_t GetCandidateBuilderConfig() const override;

    virtual target_properties_t GetTargetProperties() const override;
};

}}} // namespace ant::expconfig::setup
