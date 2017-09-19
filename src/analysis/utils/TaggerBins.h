#pragma once

#include "base/std_ext/math.h"
#include "base/interval.h"
#include "expconfig/detectors/Tagger.h"

#include "base/ParticleType.h"

#include "base/std_ext/math.h"
#include "base/std_ext/string.h"

#include "TH1D.h"

#include <iostream>
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

    static std::vector<IntervalI> EPTBinning(const unsigned nGroupedTaggerBins = 6) // this leadds to 8 taggerbins
    {
        std::vector<IntervalI> tBinning;
        for (auto i = 0 ; i < 47 ; i += nGroupedTaggerBins)
        {
            tBinning.emplace_back(i, i + nGroupedTaggerBins);
        }

        return tBinning;
    }


    static TH1D* ChannelToEg(TH1D* histCh, const tagger_t& tagger)
    {
        const std::string name = std_ext::formatter() << histCh->GetName() << "_EgBinned";
        const auto nCh = tagger->GetNChannels();

        std::vector<double> lowedges(nCh+1);
        for (auto ch = 0u; ch < nCh; ++ch)
            lowedges.at(nCh-ch-1) = tagger->GetPhotonEnergy(ch) - tagger->GetPhotonEnergyWidth(ch) / 2.;

        lowedges.at(nCh) = tagger->GetPhotonEnergy(0) + tagger->GetPhotonEnergyWidth(0) / 2.;


        auto hEg = new TH1D(name.c_str(),histCh->GetTitle(),nCh,lowedges.data());

        hEg->SetXTitle("E_{#gamma} [MeV]");
        hEg->SetYTitle(histCh->GetYaxis()->GetTitle());

        for (auto bin = 1u ;  bin <= nCh ; ++bin)
        {
            hEg->Fill(tagger->GetPhotonEnergy(nCh-bin),histCh->GetBinContent(bin));
            hEg->SetBinError(bin,histCh->GetBinError(bin));
        }
        return hEg;
    }
};



}}} // namespace ant::analysis::utils
