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

  // ReconstructHook
  virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;

  // Physics interface
  void ProcessEvent(const Event &event) override;
  void Finish() override;
  void ShowResult() override;

  // CalibrationUpdate_traits interface
  virtual std::list<TID> GetChangePoints() const override;
  virtual void Update(const TID &id) override;
};

}} // namespace ant::calibration
