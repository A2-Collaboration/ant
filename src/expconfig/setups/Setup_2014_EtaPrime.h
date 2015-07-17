
#include "ExpConfig.h"
#include "detectors/CB.h"
#include "detectors/TAPS.h"
#include "detectors/Trigger.h"
#include "unpacker/UnpackerAcqu.h"

#include "base/std_ext.h"

#include "reconstruct/Reconstruct_traits.h"
#include "calibration/modules/EnergyInvariantMass.h"
#include "calibration/modules/TimingCATCH.h"
#include "calibration/modules/IntegralSADC.h"

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EtaPrime :
    public ExpConfig::Module,
    public ExpConfig::Reconstruct,
    public UnpackerAcquConfig
{
public:
  void AddDetector(const std::shared_ptr<Detector_t>& detector) {
    detectors.push_back(detector);
  }

  template<typename T, typename... Args>
  void AddDetector(Args&&... args) {
    AddDetector(std::make_shared<T>(std::forward<Args>(args)...));
  }

  template<typename T, typename... Args>
  void AddCalibration(Args&&... args) {
    calibrations.push_back(std::make_shared<T>(std::forward<Args>(args)...));
  }

  Setup_2014_EtaPrime() {
    const auto trigger = std::make_shared<detector::Trigger>();

    AddDetector(trigger);
    AddDetector<detector::CB>();
    AddDetector<detector::TAPS_2013>(false); // no Cherenkov

    // the order of the calibrations is important
    AddCalibration<calibration::TimingCATCH>(Detector_t::Type_t::CB, trigger->Reference_CATCH_CBCrate);
    AddCalibration<calibration::IntegralSADC>(Detector_t::Type_t::CB);
  }

  virtual std::list< std::shared_ptr< CalibrationApply_traits > > GetCalibrations() const override {
    return calibrations;
  }

  virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const override {
    return detectors;
  }

  bool Matches(const THeaderInfo& header) const override {
    // check that all detectors match
    for(const auto& detector : detectors) {
      const auto& ptr = dynamic_pointer_cast<ExpConfig::Base, Detector_t>(detector);
      if(ptr == nullptr)
        continue;
      if(!ptr->Matches(header))
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
