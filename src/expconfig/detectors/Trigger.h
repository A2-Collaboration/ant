#ifndef DETECTORS_TRIGGER_H
#define DETECTORS_TRIGGER_H

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"
#include <stdexcept>

namespace ant {
namespace expconfig {
namespace detector {
struct Trigger :
    Detector_t,
    UnpackerAcquConfig
{

  Trigger() : Detector_t(Detector_t::Type_t::Trigger) {}

  virtual TVector3 GetPosition(unsigned) const override {
    throw std::runtime_error("This detector knows nothing about positions.");
  }

  virtual bool Matches(const THeaderInfo&) const override {
    return true;
  }

  const LogicalChannel_t Reference_CATCH_TaggerCrate = {Type, Channel_t::Type_t::Timing, 1000};
  const LogicalChannel_t Reference_CATCH_CBCrate = {Type, Channel_t::Type_t::Timing, 1001};

  // for UnpackerAcquConfig
  virtual void BuildMappings(
      std::vector<hit_mapping_t>& hit_mappings,
      std::vector<scaler_mapping_t>&) const override {

    hit_mappings.emplace_back(
          Reference_CATCH_TaggerCrate,
          1400
          );
    hit_mappings.emplace_back(
          Reference_CATCH_CBCrate,
          2000
          );

    /// \todo Add more data to be unpacked from trigger system
  }





};


}}} // namespace ant::expconfig::detector



#endif // DETECTORS_TRIGGER_H
