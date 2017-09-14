#pragma once

#include "base/std_ext/math.h"
#include "base/interval.h"
#include "expconfig/detectors/Tagger.h"

#include <ostream>
#include <vector>

namespace ant {
namespace analysis {
namespace utils {

struct TaggerBins {

    const std::vector<IntervalD> Bins;

    using  tagger_t = std::shared_ptr<TaggerDetector_t>;
    TaggerBins(const tagger_t& tagger, const size_t nBins):
        Bins(
            [&tagger,nBins]()
    {
        std::vector<IntervalD> bins;
        const auto nTaggCh = tagger->GetNChannels();
        const auto eMax = tagger->GetPhotonEnergy(0);
        const auto eMin = tagger->GetPhotonEnergy( nTaggCh - 1);
        const auto ebinWidth = (eMax - eMin )/ nBins;
        for (auto bin = 0u; bin < nBins; ++bin)
            bins.emplace_back(eMin + (bin * ebinWidth),eMin + ((bin+1) * ebinWidth));
        return bins;
    }()
    ){}

};



}}} // namespace ant::analysis::utils
