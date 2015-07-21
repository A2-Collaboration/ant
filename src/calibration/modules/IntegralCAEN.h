#ifndef ANT_CALIBRATION_INTEGRALCAEN_H
#define ANT_CALIBRATION_INTEGRALCAEN_H

#include "Calibration.h"
#include "expconfig/Detector_t.h"

class TH1;

namespace ant {
namespace calibration {

class IntegralCAEN : public Calibration::Module {


public:
  IntegralCAEN(
      Detector_t::Type_t detectorType
      ) :
    Calibration::Module(
      std_ext::formatter()
      << "IntegralCAEN_"
      << Detector_t::ToString(detectorType)
         ),
    DetectorType(detectorType)
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
private:
  static std::vector<double> convert(const std::vector<std::uint8_t>& rawData);
};

}}  // namespace ant::calibration

#endif
