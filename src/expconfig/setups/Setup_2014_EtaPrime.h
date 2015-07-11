
#include "ExpConfig.h"
#include "detectors/CB.h"
#include "detectors/TAPS.h"

#include "unpacker/UnpackerAcqu.h"


namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EtaPrime :
    public ExpConfig::Module,
    public UnpackerAcquConfig
{
public:
  Setup_2014_EtaPrime() {
    detectors.emplace_back(new detector::CB());
    detectors.emplace_back(new detector::TAPS());
  }

  bool Matches(const THeaderInfo& header) const override {
    return true;
  }

  void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                     std::vector<scaler_mapping_t>& scaler_mappings) const override
  {
    // this setup simply asks its underlying
    // detectors for the mappings
    for(const auto& detector : detectors) {
      const UnpackerAcquConfig* cfg
          = dynamic_cast<const UnpackerAcquConfig*>(detector.get());
      if(cfg == nullptr)
        continue;
      cfg->BuildMappings(hit_mappings, scaler_mappings);
    }
    // you may tweak the mapping at this location here
  }

private:
  std::list< std::unique_ptr<Detector_t> > detectors;
};

}}} // namespace ant::expconfig::setup
