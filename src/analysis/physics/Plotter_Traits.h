#pragma once
#include "analysis/plot/HistogramFactory.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include <functional>
#include <memory>

namespace ant {
namespace analysis {



struct Plotter_Trait {

    Plotter_Trait(WrapTFileInput& input, const HistogramFactory& HistFactory, OptionsPtr otps) {}

    virtual long long GetNumEntries() const =0;
    virtual bool ProcessEntry(const long long entry) =0;
    virtual void Finish() {}
    virtual void ShowResult() {}

    virtual ~Plotter_Trait() {}

};

using Plotter_Maker_t = std::function< std::unique_ptr<Plotter_Trait>(WrapTFileInput& input, const HistogramFactory& HistFactory, OptionsPtr otps) >;

}
}
