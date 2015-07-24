#pragma once

#include "Calibration.h"

class TH1;

namespace ant {
namespace calibration {

class EnergyInvariantMass : public Calibration::Module {
protected:
  TH1* ggIM = nullptr;

public:
  EnergyInvariantMass();

  // CalibrationApply_traits interface
  virtual void ApplyTo(const readhits_t& hits) override;
  virtual void ApplyTo(event_ptr&) override {}

  // Physics interface
  void ProcessEvent(const Event &event) override;
  void Finish() override;
  void ShowResult() override;

  // CalibrationUpdate_traits interface
  void BuildRanges(std::list<TID> &ranges) override;
  void Update(const TID &id) override;
};

}} // namespace ant::calibration
