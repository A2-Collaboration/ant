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

class Setup_2010_03_Base : public Setup
{
protected:
    const bool MCTaggerHits;
    const bool cherenkovInstalled;
    const std::shared_ptr<detector::Trigger_2007> Trigger;
    const std::shared_ptr<detector::Tagger_2010_03> Tagger;
    const std::shared_ptr<detector::CB> CB;
    const std::shared_ptr<detector::PID_2009_07> PID;
    const std::shared_ptr<detector::TAPS_2009_03> TAPS;
    const std::shared_ptr<detector::TAPSVeto_2009_03> TAPSVeto;
public:

    Setup_2010_03_Base(const std::string& name, OptionsPtr opt);

    virtual double GetElectronBeamEnergy() const override;

    virtual candidatebuilder_config_t GetCandidateBuilderConfig() const override;
};

}}} // namespace ant::expconfig::setup
