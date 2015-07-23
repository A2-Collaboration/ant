#pragma once

#include "Calibration.h"
#include "expconfig/Detector_t.h"

class TH1;

namespace ant {
namespace calibration {

class TimingTAPS : public Calibration::Module {


public:
  TimingTAPS() :
    Calibration::Module("TimingTAPS")
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

private:
  static std::vector<double> convert(const std::vector<std::uint8_t>& rawData);
};

}}  // namespace ant::calibration
