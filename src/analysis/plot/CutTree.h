#pragma once

#include "analysis/plot/HistogramFactories.h"
#include "root-addons/analysis_codes/hstack.h"
#include "base/Tree.h"

#include "TH1.h"

#include <vector>
#include <list>
#include <map>
#include <string>
#include <functional>

namespace ant {
namespace analysis {
namespace plot {

/**
 * @brief Please see progs/EtapOmegaG_plot.cc how to use cuttree
 */
namespace cuttree {

template<typename Fill_t>
struct Cut_t {
    using Passes_t = std::function<bool(const Fill_t&)>;
    Cut_t(const std::string& name,
          Passes_t passes = [] (const Fill_t&) { return true; }) : Name(name), Passes(passes) {}
    const std::string Name;
    const Passes_t Passes;
};

template<typename Fill_t>
using MultiCut_t = std::vector<Cut_t<Fill_t>>;

template<typename Fill_t>
using Cuts_t = std::list<MultiCut_t<Fill_t>>;

template<typename ToFill_t, typename FromFill_t>
Cuts_t<ToFill_t> ConvertCuts(const Cuts_t<FromFill_t>& from_cuts) {
    Cuts_t<ToFill_t> to_cuts;
    for(const auto& from_multicut : from_cuts) {
        to_cuts.emplace_back();
        auto& to_multicut = to_cuts.back();
        for(const auto& from_cut : from_multicut) {
            auto passes = from_cut.Passes;
            to_multicut.emplace_back(from_cut.Name, [passes] (const ToFill_t& f) { return passes(f); });
        }
    }
    return to_cuts;
}

template<typename Hist_t>
struct Node_t {
    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    HistogramFactory HistFac;
    Hist_t Hist;
    typename Cut_t<Fill_t>::Passes_t PassesCut;

    Node_t(const HistogramFactory& parentHistFac, Cut_t<Fill_t> cut) :
        HistFac(cut.Name, parentHistFac, cut.Name),
        Hist(HistFac),
        PassesCut(cut.Passes)
    {}

    bool Fill(const Fill_t& f) {
       if(PassesCut(f)) {
           Hist.Fill(f);
           return true;
       }
       return false;
    }
};

template<typename Hist_t>
using CutsIterator_t = typename Cuts_t<typename Hist_t::Fill_t>::const_iterator;

template<typename Hist_t>
using Tree_t = typename Tree<Node_t<Hist_t>>::node_t;

template<typename Hist_t>
void Build(Tree_t<Hist_t> cuttree, CutsIterator_t<Hist_t> first, CutsIterator_t<Hist_t> last) {
    if(first == last)
        return;
    const auto& multicut = *first;
    for(const auto& cut : multicut) {
        auto daughter = cuttree->CreateDaughter(cuttree->Get().HistFac, cut);
        Build<Hist_t>(daughter, std::next(first), last);
    }
}

template<typename Hist_t, typename Fill_t = typename Hist_t::Fill_t>
Tree_t<Hist_t> Make(HistogramFactory histFac, const std::string& name, const Cuts_t<Fill_t>& cuts) {
    auto cuttree = Tree<Node_t<Hist_t>>::MakeNode(histFac, Cut_t<Fill_t>{name});
    Build<Hist_t>(cuttree, cuts.begin(), cuts.end());
    return cuttree;
}

template<typename Hist_t, typename Fill_t = typename Hist_t::Fill_t>
void Fill(Tree_t<Hist_t> cuttree, const Fill_t& f) {
    if(cuttree->Get().Fill(f)) {
        for(const auto& d : cuttree->Daughters()) {
            Fill<Hist_t>(d, f);
        }
    }
}

// returns string to be used as draw option
struct HistMod_t : std::function<std::string(TH1*)> {
    // use constructors
    using std::function<std::string(TH1*)>::function;

    HistMod_t() : std::function<std::string(TH1*)>([] (TH1*) { return ""; }) {}

    static Color_t GetColor(unsigned i) {
        const std::vector<Color_t> colors = {kGreen+1, kBlue, kYellow+1, kMagenta, kCyan, kOrange, kSpring+10,};
        return colors[i % colors.size()];
    }

    static HistMod_t MakeLine(const Color_t color, short linewidth = 1) {
        return [color, linewidth] (TH1* h) {
            h->SetLineColor(color);
            h->SetLineWidth(linewidth);
            h->SetMarkerSize(1);
            h->SetMarkerColor(color);
            return "";
        };
    }

    static HistMod_t MakeDataPoints(const Color_t color, short linewidth = 1) {
        return [color, linewidth] (TH1* h) {
            h->SetLineColor(color);
            h->SetLineWidth(linewidth);
            h->SetMarkerColor(color);
            h->SetMarkerStyle(kDot);
            return "E"; // draw error bars
        };

    }

};

template<typename Hist_t>
struct StackedHists_t {
protected:
    StackedHists_t(const HistogramFactory& histFac) :
        HistFac(histFac),
        H("h", HistFac)
    {
    }


    const Hist_t& GetHist(unsigned key,
                          const std::string& name = "",
                          HistMod_t histmod = {}) {
        auto it_hist = Hists.lower_bound(key);
        if(it_hist == Hists.end() || it_hist->first != key) {
            if(name.empty())
                throw std::runtime_error("Cannot create hist with empty name");
            it_hist = Hists.emplace_hint(it_hist, key, WrappedHist_t{
                                             name,
                                             H,
                                             histmod
                                         });

            AddToStack(it_hist->second);
        }
        return it_hist->second.Hist;
    }



private:

    HistogramFactory HistFac; // directory for stacks and H itself
    HistogramFactory H; // directory for mctrue splitted histfacs

    struct WrappedHist_t {
        Hist_t    Hist;
        HistMod_t Modify;
        WrappedHist_t(const std::string& name,
                      const HistogramFactory& parentHistFac,
                      HistMod_t mod) :
            Hist(HistogramFactory(name, parentHistFac, name)), Modify(mod) {}
    };

    std::map<unsigned, WrappedHist_t> Hists;
    std::vector<ant::hstack*> Stacks;

    void AddToStack(const WrappedHist_t& hist) {
        const auto histptrs = hist.Hist.GetHists();

        if(Stacks.empty()) {
            for(auto h : histptrs) {
                const std::string name = h->GetName();
                Stacks.emplace_back(
                            HistFac.make<ant::hstack>(
                                // use the parent histFac here!
                                name,
                                H.GetTitlePrefix()+": "+name,
                                ant::hstack::options_t::all_enabled
                                ));
            }
        }
        assert(histptrs.size() == Stacks.size());

        for(size_t i=0;i<Stacks.size();i++) {
            TH1* h = histptrs[i];
            *Stacks[i] << drawoption(hist.Modify(h)) << h;
        }
    }

};

} // namespace cuttree

}}}
