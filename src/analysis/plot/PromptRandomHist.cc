#include "PromptRandomHist.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::PromptRandom;
using namespace std;

void Hist1::MakeHistograms(const HistogramFactory& factory, const string& name, const string& title, const BinSettings& bins, const string& xtitle, const string& ytitle) {

    HistogramFactory myFactory(name, factory, title);

    prompt     = myFactory.makeTH1D("Prompt",    xtitle, ytitle, bins, "prompt");
    random     = myFactory.makeTH1D("Random",    xtitle, ytitle, bins, "random");
    subtracted = myFactory.makeTH1D("Subtracted",xtitle, ytitle, bins, "subtracted");

    subtracted->Sumw2();
}

void Hist2::MakeHistograms(const HistogramFactory& factory, const string& name, const string& title, const BinSettings& xbins, const BinSettings& ybins, const string& xtitle, const string& ytitle) {

    HistogramFactory myFactory(name, factory, title);

    prompt     = myFactory.makeTH2D("Prompt",    xtitle, ytitle, xbins, ybins, "prompt");
    random     = myFactory.makeTH2D("Random",    xtitle, ytitle, xbins, ybins, "random");
    subtracted = myFactory.makeTH2D("Subtracted",xtitle, ytitle, xbins, ybins, "subtracted");

    subtracted->Sumw2();
}

void Switch::update_ratio() {
    double p = promptw.Area();
    double r = randomw.Area();
    if(r <= 0.0) {
        ratio = 1.0;
    } else {
        ratio = p/r;
    }
}

Switch::Switch(std::initializer_list<Switch::interval_t> p, std::initializer_list<Switch::interval_t> r):
    promptw(p),
    randomw(r) { update_ratio(); }

Switch::Switch(Switch::windows_t p, Switch::windows_t r):
    promptw(p),
    randomw(r) { update_ratio(); }

void Switch::AddPromptRange(const Switch::interval_t& i) {
    promptw.emplace_back(i);
    update_ratio();
}

void Switch::AddRandomRange(const Switch::interval_t& i) {
    randomw.emplace_back(i);
    update_ratio();
}

void Switch::SetTaggerHit(const T tagtime) {

    if(randomw.Contains(tagtime)) {
        rpcase = Case::Random;
        fillw = -Ratio();
    } else if(promptw.Contains(tagtime)) {
        rpcase = Case::Prompt;
        fillw = 1.0;
    } else {
        rpcase = Case::Outside;
        fillw = 0.0;
    }

}
