#pragma once

#include "expconfig/Detector_t.h"

#include <memory>
#include <map>
#include <list>

namespace ant {

class TID;
class TDetectorRead;
class TDetectorReadHit;
class TEvent;



/**
 * @brief The CalibrationApply_traits class
 * Applies calibration factors to detector reads
 */
class CalibrationApply_traits {
public:
    using readhits_t = std::map< Detector_t::Type_t, std::list< TDetectorReadHit* > >;
    using event_ptr = std::unique_ptr<TEvent>;

    virtual void ApplyTo(TDetectorRead& detectorRead, const readhits_t& hits) = 0;
    virtual void ApplyTo(event_ptr& event) = 0;
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
