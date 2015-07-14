#ifndef ANT_RECONSTRUCT_H
#define ANT_RECONSTRUCT_H

#include <memory>

namespace ant {

class THeaderInfo;
class TDetectorRead;
class TEvent;

class Reconstruct {
public:
  // You can only create the reconstruct machinery
  // if it's able to find its config. For this, it needs the
  // the headerInfo
  Reconstruct(const THeaderInfo& headerInfo);

  // this method converts a TDetectorRead
  // into a calibrated TEvent
  std::unique_ptr<TEvent> DoReconstruct(std::unique_ptr<TDetectorRead>& read);
};

}

#endif // ANT_RECONSTRUCT_H
