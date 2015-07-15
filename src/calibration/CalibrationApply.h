#ifndef CALIBRATIONAPPLY_H
#define CALIBRATIONAPPLY_H

#include <memory>
#include <list>

#include "unpacker/tree/TDataRecord.h"
#include "base/interval.h"

namespace ant {

class TDetectorRead;

/**
 * @brief The CalibrationApply_traits class
 * Applies calibration factors to detector reads
 */
class CalibrationApply_traits {
public:
  virtual void ApplyTo(std::unique_ptr<TDetectorRead>& detectorRead) = 0;
};


/**
 * @brief The CalibrationUpdate_traits class
 *
 */
class CalibrationUpdate_traits {
  virtual void BuildRanges(std::list<TDataRecord::ID_t>& ranges) = 0;
  virtual void Update(const TDataRecord::ID_t& id) = 0;
};

} // namespace ant

#endif // CALIBRATIONAPPLY_H
