
#include "ExpConfig.h"
#include "detectors/CB.h"
#include "detectors/TAPS.h"
#include "unpacker/UnpackerAcqu.h"

#include "base/std_ext.h"


namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EtaPrime :
    public ExpConfig::Module,
    public UnpackerAcquConfig
{
public:
  Setup_2014_EtaPrime() {
    detectors.push_back(std_ext::make_unique<detector::CB>());
    detectors.push_back(std_ext::make_unique<detector::TAPS_2013>(false)); // no Cherenkov
  }

  virtual void GetCalibrations() const override {}

  bool Matches(const THeaderInfo& header) const override {
    // check that all detectors match
    for(const auto& detector : detectors) {
      const ExpConfig::Base* cfg
          = dynamic_cast<const ExpConfig::Base*>(detector.get());
      if(cfg == nullptr)
        continue;
      if(!cfg->Matches(header))
        return false;
    }
    /// \todo Make beamtime match stricter than just detectors
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
    // for example, ignore elements
  }

private:
  std::list< std::unique_ptr<Detector_t> > detectors;
};

}}} // namespace ant::expconfig::setup
