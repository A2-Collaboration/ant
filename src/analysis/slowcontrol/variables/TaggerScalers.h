#pragma once

#include "Variable.h"

#include "base/Detector_t.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace variable {

struct TaggerScalers : Variable {

    virtual void Init() override;
    virtual std::list<ProcessorPtr> GetNeededProcessors() const override;

    /**
     * @brief Get returns the tagger current scaler frequencies
     * @return
     */
    std::vector<double> Get() const;

    /**
     * @brief GetCounts returns the counts for current scalar block
     * @return
     */
    std::vector<int64_t> GetCounts() const;

    double GetTaggerOr() const;

    double GetTaggerOrAsSum() const;

protected:

    enum class mode_t {
        EPT_2014,
        Tagger
    };

    mode_t mode;
    unsigned nChannels;

};

}}}} // namespace ant::analysis::slowcontrol::processor
