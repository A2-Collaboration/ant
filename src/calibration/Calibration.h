#ifndef ANT_CALIBRATION_H
#define ANT_CALIBRATION_H

#include "analysis/physics/Physics.h"
#include "reconstruct/Reconstruct_traits.h"

namespace ant {


class Calibration {
public:

  /**
   * @brief The Module class
   * A calibration module is two things:
   * * A physics class to make histograms and determine calibration factors
   * * and a CalibrationApply class that can apply those factors to data
   * * furthermore, it may update its constants during the analysis
   */
  class Module :
      public Physics,
      public CalibrationApply_traits,
      public Updateable_traits
  {
  public:
    Module(const std::string& name): Physics(name) {}
    virtual ~Module() = default;
  };

};

} // namespace ant

#endif // ANT_CALIBRATION_H
