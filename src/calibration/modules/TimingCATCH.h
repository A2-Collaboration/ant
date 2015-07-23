#pragma once

#include "Calibration.h"
#include "expconfig/Detector_t.h"

class TH1;

namespace ant {
namespace calibration {

class TimingCATCH : public Calibration::Module {


public:
  TimingCATCH(
      Detector_t::Type_t detectorType,
      const LogicalChannel_t& referenceChannel
      ) :
    Calibration::Module(
      std_ext::formatter()
      << "TimingCATCH_"
      << Detector_t::ToString(detectorType)
         ),
    DetectorType(detectorType),
    ReferenceChannel(referenceChannel)
  {}

  // CalibrationApply_traits interface
  virtual void ApplyTo(const std::map< Detector_t::Type_t, std::list< TDetectorReadHit* > >& hits) override;
  virtual void ApplyTo(std::unique_ptr<TEvent>&) override {}


  // Physics interface
  void ProcessEvent(const Event &event) override;
  void Finish() override;
  void ShowResult() override;

  // CalibrationUpdate_traits interface
  void BuildRanges(std::list<TID>&) override {}
  void Update(const TID&) override {}

protected:
  const Detector_t::Type_t DetectorType;
  const LogicalChannel_t   ReferenceChannel;

private:
  static std::vector<double> convert(const std::vector<std::uint8_t>& rawData);
};

}}  // namespace ant::calibration
