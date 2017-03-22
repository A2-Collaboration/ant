#pragma once

#include <vector>

#include <memory>
#include <base/Interpolator.h>

#include <ostream>

class TH2D;

//---------------- copy & paste from interpolated pulls ----------------------------

namespace ant {

struct ClippedInterpolatorWrapper {
    using interpolator_ptr_t = std::unique_ptr<const ant::Interpolator2D>;

    interpolator_ptr_t interp;

    struct boundsCheck_t {
        ant::interval<double> range;
        mutable unsigned underflow = 0;
        mutable unsigned unclipped = 0;
        mutable unsigned overflow  = 0;

        double clip(double v) const;

        boundsCheck_t(const ant::interval<double> r): range(r) {}
    };
    friend std::ostream& operator<<(std::ostream& stream, const boundsCheck_t& o);

    boundsCheck_t xrange;
    boundsCheck_t yrange;

    ClippedInterpolatorWrapper(interpolator_ptr_t i);
    ClippedInterpolatorWrapper();
    ~ClippedInterpolatorWrapper();
    double GetPoint(double x, double y) const;

    void setInterpolator(interpolator_ptr_t i);

    friend std::ostream& operator<<(std::ostream& stream, const ClippedInterpolatorWrapper& o);

    static std::unique_ptr<const Interpolator2D> makeInterpolator(TH2D* hist);

    TH2D* getCheckHistogram(const int xBins, const int yBins);

};




} // namespace ant

