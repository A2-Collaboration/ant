#pragma once

#include "Processor.h"

#include "tree/TEventData.h"

#include <queue>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct TaggEff : Processor {

    virtual return_t ProcessEventData(const TEventData& recon,  physics::manager_t& manager) override;

    virtual void PopQueue() override;

    using value_t = std::vector<TaggerDetector_t::taggeff_t>;
    value_t Get() const;

private:
    std::queue<value_t> queue;
};

}}}} // namespace ant::analysis::slowcontrol::processor