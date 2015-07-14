#ifndef CALIBRATIONAPPLY_H
#define CALIBRATIONAPPLY_H

#include <memory>

namespace ant {

class TDetectorRead;

/**
 * @brief The CalibrationApply class
 * Applies calibration factors to data
 * @todo think about this and invent something nice...
 */
class CalibrationApply_traits {
public:
  virtual void ApplyTo(std::unique_ptr<TDetectorRead>& detectorRead) = 0;
};

} // namespace ant

#endif // CALIBRATIONAPPLY_H
