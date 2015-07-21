#ifndef ANT_RECONSTRUCT_TRAITS_H
#define ANT_RECONSTRUCT_TRAITS_H

#include <memory>
#include <map>

#include "expconfig/Detector_t.h"

namespace ant {

class TID;
class TDetectorReadHit;
class TEvent;



/**
 * @brief The CalibrationApply_traits class
 * Applies calibration factors to detector reads
 */
class CalibrationApply_traits {
public:
    virtual void ApplyTo(const std::map< Detector_t::Type_t, std::list< TDetectorReadHit* > >& hits) = 0;
    virtual void ApplyTo(std::unique_ptr<TEvent>& event) = 0;
};


/**
 * @brief The Updateable_traits class
 *
 */
class Updateable_traits {
public:
    virtual void BuildRanges(std::list<TID>& ranges) = 0;
    virtual void Update(const TID& id) = 0;
};

} // namespace ant

#endif // ANT_RECONSTRUCT_TRAITS_H
