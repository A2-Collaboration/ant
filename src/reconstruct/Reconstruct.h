#ifndef ANT_RECONSTRUCT_H
#define ANT_RECONSTRUCT_H

#include <memory>
#include <list>

#include "Reconstruct_traits.h"
#include "tree/TEvent.h"

namespace ant {

class THeaderInfo;
class TDetectorRead;

class Reconstruct {
public:
  // You can only create the reconstruct machinery
  // if it's able to find its config. For this, it needs the
  // some THeaderInfo object
  Reconstruct(const THeaderInfo& headerInfo);

  // this method converts a TDetectorRead
  // into a calibrated TEvent
  std::unique_ptr<TEvent> DoReconstruct(TDetectorRead& detectorRead);

private:
  std::list< std::shared_ptr<CalibrationApply_traits> >  calibrations;
  std::list< std::shared_ptr<Updateable_traits> >  updateables;


};

}

#endif // ANT_RECONSTRUCT_H
