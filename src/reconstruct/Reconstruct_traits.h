#ifndef ANT_RECONSTRUCT_TRAITS_H
#define ANT_RECONSTRUCT_TRAITS_H

#include <memory>
#include <list>

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"
#include "base/interval.h"

namespace ant {

/**
 * @brief The CalibrationApply_traits class
 * Applies calibration factors to detector reads
 */
class CalibrationApply_traits {
public:
  virtual void ApplyTo(std::unique_ptr<TDetectorRead>& detectorRead) = 0;
  virtual void ApplyTo(std::unique_ptr<TEvent>& event) = 0;
};


/**
 * @brief The Updateable_traits class
 *
 */
class Updateable_traits {
  virtual void BuildRanges(std::list<TDataRecord::ID_t>& ranges) = 0;
  virtual void Update(const TDataRecord::ID_t& id) = 0;
};

} // namespace ant

#endif // ANT_RECONSTRUCT_TRAITS_H
