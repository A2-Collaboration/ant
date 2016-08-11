#pragma once

#include "Variable.h"

#include "base/Detector_t.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace variable {

struct PhotonFlux : Variable {

    virtual std::list<ProcessorPtr> GetNeededProcessors() override;

    double Get() const;

protected:

    enum class mode_t {
        PbGlass
    };

    mode_t mode;

};

}}}} // namespace ant::analysis::slowcontrol::processor
