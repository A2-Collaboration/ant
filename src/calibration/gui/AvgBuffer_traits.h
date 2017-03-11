#pragma once

#include "base/interval.h"
#include "tree/TID.h"

#include <memory>

namespace ant {
namespace calibration {
namespace gui {

template<typename Hist>
class AvgBuffer_traits {
public:
    virtual void Peek(const interval<TID>& range) {
        total_length += range.Stop().Lower - range.Start().Lower;
        ++total_n;
    }
    virtual void Push(std::shared_ptr<Hist> hist, const interval<TID>& range) =0;
    virtual bool Empty() const =0;
    virtual void Flush() =0;
    virtual void Next() =0;

    virtual const Hist& CurrentHist() const =0;
    virtual const interval<TID>& CurrentRange() const =0;

    virtual ~AvgBuffer_traits() = default;
protected:
    double   total_length = 0;
    unsigned total_n = 0;
};

}}} // namespace ant::calibration::gui