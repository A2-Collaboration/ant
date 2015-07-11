
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
  bool Matches(const THeaderInfo& header) const override {
    return true;
  }

  void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                     std::vector<scaler_mapping_t>& scaler_mappings) override
  {
    hit_mapping_t hit_map;
    hit_map.LogicalChannel = {Detector_t::Type_t::Trigger, Channel_t::Type_t::Counter, 31};
    hit_map.RawChannels.push_back(400);
//    map.LogicalElement = {Detector_t::Trigger, ChannelType_t::Counter, 1};
//    map.RawChannels.push_back(1853);
    hit_mappings.push_back(hit_map);

    scaler_mapping_t scaler_map;
    scaler_map.LogicalChannel = {Detector_t::Type_t::Trigger, Channel_t::Type_t::Scaler, 17};
    scaler_map.RawChannels.push_back(0);
    scaler_mappings.push_back(scaler_map);
    scaler_map.SlowControlName = "Some trigger scaler";
    scaler_mappings.push_back(scaler_map);
  }


};

}}} // namespace ant::expconfig::setup
