#include "Setup.h"

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Common base class for all Setups of the 2014 Ent Point Tagger beam times
 */
class Setup_2007_Base : public Setup
{
    const bool MCTaggerHits;

public:

    Setup_2007_Base(const std::string& name, OptionsPtr opt);

    virtual double GetElectronBeamEnergy() const override;

    virtual ExpConfig::Setup::candidatebuilder_config_t GetCandidateBuilderConfig() const override;
};

}}} // namespace ant::expconfig::setup
