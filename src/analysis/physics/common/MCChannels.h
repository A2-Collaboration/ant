#include "analysis/physics/Physics.h"

#include "base/ParticleType.h"

#include <map>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief Phyics class for counting occurance of photoproduction channels in MC Data.
 */
class MCChannels : public Physics {

    using Counter_t = std::map<std::string, unsigned>;

    Counter_t counter_production;
    unsigned total  = 0;
    unsigned noTree = 0;

    TH1D* h_production = nullptr;

    TH1D* h_database = nullptr;
    TH2D* h_database_taggch = nullptr;

public:
    MCChannels(const std::string& name, OptionsPtr opts);
    virtual ~MCChannels();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;

};

}}}
