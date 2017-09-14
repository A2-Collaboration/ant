#pragma once

#include "base/std_ext/math.h"
#include "base/interval.h"
#include "expconfig/detectors/Tagger.h"

#include "base/ParticleType.h"

#include "base/std_ext/math.h"

#include <ostream>
#include <vector>

namespace ant {
namespace analysis {
namespace utils {

struct TaggerBins {
    using tagger_t = std::shared_ptr<TaggerDetector_t>;
    using ranges_t = std::vector<IntervalD>;

    static ranges_t MakeEgBins(const tagger_t& tagger, const size_t nBins)
    {
        ranges_t bins;
        const auto nTaggCh = tagger->GetNChannels();
        const auto eMax = tagger->GetPhotonEnergy(0);
        const auto eMin = tagger->GetPhotonEnergy( nTaggCh - 1);
        const auto ebinWidth = (eMax - eMin )/ nBins;
        for (auto bin = 0u; bin < nBins; ++bin)
            bins.emplace_back(eMin + (bin * ebinWidth),eMin + ((bin+1) * ebinWidth));
        return bins;
    }

    static ranges_t MakeWBins(const tagger_t& tagger, const size_t nBins)
    {
        auto mevBins = MakeEgBins(tagger,nBins);
        ranges_t WBins;
        std::transform(mevBins.begin(),mevBins.end(),std::back_inserter(WBins),
                      [](IntervalD& range)
        {
            auto W = [](const double eg)
            {
                return sqrt(std_ext::sqr(eg + ParticleTypeDatabase::Proton.Mass())
                            - std_ext::sqr(eg));
            };

            return IntervalD(W(range.Start()),W(range.Stop()));
        });
        return WBins;
    }
};



}}} // namespace ant::analysis::utils
