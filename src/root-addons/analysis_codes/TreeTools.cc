#include "TreeTools.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCut.h"
#include "TDirectory.h"

using namespace std;

TH1*Draw(TTree* tree, const string& formula, const TCut& cut, const int bins, const double min, const double max) {
    static unsigned n = 0;
    const char* hname = Form("h1d_%d", n++);

    tree->Draw(Form("%s>>%s(%d,%lf,%lf)",formula.c_str(),hname,bins,min,max),cut);
    TH1* h = NULL;
    gDirectory->GetObject(hname, h);
    return h;
}

TH2*Draw(TTree* tree, const string& formula, const TCut& cut, const int xbins, const double xmin, const double xmax, const int ybins, const double ymin, const double ymax) {
    static unsigned n = 0;
    const char* hname = Form("h2d_%d", n++);

    tree->Draw(Form("%s>>%s(%d,%lf,%lf,%d,%lf,%lf)",formula.c_str(),hname,xbins,xmin,xmax,ybins,ymin,ymax),cut,"colz");
    TH2* h = NULL;
    gDirectory->GetObject(hname, h);
    return h;
}

TH3* Draw(TTree* tree, const string& formula, const TCut& cut, const ant::analysis::BinSettings& xbins, const ant::analysis::BinSettings& ybins, const ant::analysis::BinSettings& zbins)
{
    static unsigned n = 0;
    const char* hname = Form("h3d_%d", n++);

    tree->Draw(
                Form(
                    "%s>>%s(%d,%lf,%lf,%d,%lf,%lf,%d,%lf,%lf)",
                    formula.c_str(),
                    hname,
                    xbins.Bins(), xbins.Start(), xbins.Stop(),
                    ybins.Bins(), ybins.Start(), ybins.Stop(),
                    zbins.Bins(), zbins.Start(), zbins.Stop()
                    )
                ,cut);
    TH3* h = NULL;
    gDirectory->GetObject(hname, h);
    return h;
}

TH1*Draw(TTree* tree, const string& formula, const TCut& cut, const string& xtitle, const string& ytitle, const ant::analysis::BinSettings& xbins)
{
    auto h = Draw(tree, formula, cut, int(xbins.Bins()), xbins.Start(), xbins.Stop());
    h->SetXTitle(xtitle.c_str());
    h->SetYTitle(ytitle.c_str());

    return h;
}

TH2*Draw(TTree* tree, const string& formula, const TCut& cut, const string& xtitle, const string& ytitle, const ant::analysis::BinSettings& xbins, const ant::analysis::BinSettings& ybins)
{
    auto h = Draw(tree, formula, cut, int(xbins.Bins()), xbins.Start(), xbins.Stop(), int(ybins.Bins()), ybins.Start(), ybins.Stop());
    h->SetXTitle(xtitle.c_str());
    h->SetYTitle(ytitle.c_str());

    return h;
}
