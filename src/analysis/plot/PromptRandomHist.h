#pragma once

#include "base/interval.h"
#include "base/piecewise_interval.h"
#include "HistogramFactory.h"

#include "TH1D.h"
#include "TH2D.h"

#include <string>

namespace ant {
namespace analysis {
namespace PromptRandom {

enum class Case {
    Prompt,
    Random,
    Outside
};


class Switch {
public:
    using T = double;
    using windows_t  = ant::PiecewiseInterval<T>;
    using interval_t = windows_t::interval_t;

protected:
    windows_t promptw = {};
    windows_t randomw = {};
    double ratio = 1.0;  // == all prompt

    void update_ratio();

    Case rpcase = Case::Prompt;
    double fillw = 1.0;

public:
    Switch() {}
    Switch(std::initializer_list<interval_t> p, std::initializer_list<interval_t> r);
    Switch(windows_t p, windows_t r);

    // Since Switch and Hist1/2 are loosely bound
    // we do not allow moves/copies after Switch was created
    Switch(const Switch&) = delete;
    Switch& operator=(Switch&) = delete;
    Switch(Switch&&) = delete;
    Switch& operator=(Switch&&) = delete;



    /**
     * @brief
     * @return
     */
    double Ratio() const { return ratio; }

    double FillWeight() const { return fillw; }

    Case State() const { return rpcase; }

    void AddPromptRange(const interval_t& i);

    void AddRandomRange(const interval_t& i);

    void SetTaggerHit(const T tagtime);


};


struct Hist1 {
    TH1D* prompt;
    TH1D* random;
    TH1D* subtracted;

    const Switch& d;

    Hist1(Switch& D):
        d(D) {}

    void MakeHistograms(const HistogramFactory& factory,
                        const std::string& name,
                        const std::string& title,
                        const BinSettings& bins,
                        const std::string& xtitle,
                        const std::string& ytitle);

    void Fill(const double x) {

        subtracted->Fill(x, d.FillWeight());

        if(d.State() == Case::Prompt) {
            prompt->Fill(x);
        } else if(d.State() == Case::Random) {
            random->Fill(x);
        }

    }
};

struct Hist2 {
    TH2D* prompt;
    TH2D* random;
    TH2D* subtracted;

    const Switch& d;

    Hist2(Switch& D):
        d(D) {}

    void MakeHistograms(const HistogramFactory& factory,
                        const std::string& name,
                        const std::string& title,
                        const BinSettings& xbins,
                        const BinSettings& ybins,
                        const std::string& xtitle,
                        const std::string& ytitle);

    void Fill(const double x, const double y) {

        subtracted->Fill(x, y, d.FillWeight());

        if(d.State() == Case::Prompt) {
            prompt->Fill(x, y);
        } else if(d.State() == Case::Random) {
            random->Fill(x, y);
        }

    }
};

}
}
}
