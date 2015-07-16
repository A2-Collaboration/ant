
#include "ExpConfig.h"
#include "detectors/CB.h"
#include "detectors/TAPS.h"
#include "unpacker/UnpackerAcqu.h"

#include "base/std_ext.h"

#include "calibration/CalibrationApply.h"
#include "calibration/TestCalCB.h"

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EtaPrime :
    public ExpConfig::Module,
    public ExpConfig::Reconstruct,
    public UnpackerAcquConfig
{
public:
  Setup_2014_EtaPrime() {
    detectors.push_back(std::make_shared<detector::CB>());
    detectors.push_back(std::make_shared<detector::TAPS_2013>(false)); // no Cherenkov
    calibrations.push_back(std::make_shared<calibration::TestCalCB>());
  }

  virtual std::list< std::shared_ptr< CalibrationApply_traits > > GetCalibrations() const override {
    return calibrations;
  }

  virtual std::list< std::shared_ptr< CalibrationUpdate_traits > > GetUpdateables() const override {
    return {};
  }

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
                     std::vector<scaler_mapping_t>& scaler_mappings) const
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
  std::list< std::shared_ptr<Detector_t> > detectors;
  std::list< std::shared_ptr<CalibrationApply_traits> > calibrations;
};

}}} // namespace ant::expconfig::setup
