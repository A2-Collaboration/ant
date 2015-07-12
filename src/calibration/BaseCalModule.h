#ifndef BASECALMODULE_H
#define BASECALMODULE_H

#include "analysis/physics/Physics.h"

namespace ant {
namespace calibration {


/**
 * @brief The CalibrationApply class
 * Applies calibration factors to data
 * @todo think about this and invent something nice...
 */
class CalibrationApply_traits {
public:
    virtual void ApplyTo() = 0;
};

/**
 * @brief The BaseCalibrationModule class
 * A Calibration module is two things:
 * * A physics class to make histograms and determine calibration factors
 * * and a CalibrationApply class that can apply those factors to data
 */
class BaseCalibrationModule: public ant::Physics, public CalibrationApply_traits {
public:
    BaseCalibrationModule(const std::string& name): Physics(name) {}
    virtual ~BaseCalibrationModule() = default;
};

}
}

#endif
