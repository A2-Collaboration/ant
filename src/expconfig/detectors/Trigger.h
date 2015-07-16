#ifndef DETECTORS_TRIGGER_H
#define DETECTORS_TRIGGER_H

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct Trigger :
    Detector_t,
    UnpackerAcquConfig
{

  Trigger() : Detector_t(Detector_t::Type_t::Trigger) {}

  virtual bool Matches(const THeaderInfo&) const override {
    return true;
  }

  // for UnpackerAcquConfig
  virtual void BuildMappings(
      std::vector<hit_mapping_t>& hit_mappings,
      std::vector<scaler_mapping_t>&) const override {

    hit_mapping_t refCATCH1;
    refCATCH1.LogicalChannel = {Type, Channel_t::Type_t::Timing, 1000};
    refCATCH1.RawChannels.push_back(1400);
    hit_mappings.emplace_back(std::move(refCATCH1));

    hit_mapping_t refCATCH2;
    refCATCH2.LogicalChannel = {Type, Channel_t::Type_t::Timing, 1001};
    refCATCH2.RawChannels.push_back(2000);
    hit_mappings.emplace_back(std::move(refCATCH2));

    /// \todo Add more data to be unpacked from trigger system
  }





};


}}} // namespace ant::expconfig::detector



#endif // DETECTORS_TRIGGER_H
