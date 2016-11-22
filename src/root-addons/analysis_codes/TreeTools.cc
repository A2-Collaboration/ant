#include "TreeTools.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCut.h"
#include "TDirectory.h"
#include "base/std_ext/string.h"

using namespace std;
using namespace ant::std_ext;
using namespace ant::analysis;


TH1*Draw(TTree* tree, const string& formula, const TCut& cut, const string& xtitle, const string& ytitle, const ant::BinSettings& xbins, const std::string& name)
{
    TH1* h = new TH1D(name.c_str(),"",int(xbins.Bins()), xbins.Start(), xbins.Stop());

    tree->Draw(Form("%s>>%s",formula.c_str(),name.c_str()),cut);

    h->SetXTitle(xtitle.c_str());
    h->SetYTitle(ytitle.c_str());

    return h;
}

TH2*Draw(TTree* tree, const string& formula, const TCut& cut, const string& xtitle, const string& ytitle, const ant::BinSettings& xbins, const ant::BinSettings& ybins, const std::string& name)
{
    TH2* h = new TH2D(name.c_str(),"",int(xbins.Bins()), xbins.Start(), xbins.Stop(), int(ybins.Bins()), ybins.Start(), ybins.Stop());
    tree->Draw(Form("%s>>%s",formula.c_str(), name.c_str()), cut,"colz");
    h->SetXTitle(xtitle.c_str());
    h->SetYTitle(ytitle.c_str());

    return h;
}

TH3* Draw(TTree* tree, const string& formula, const TCut& cut, const ant::BinSettings& xbins, const ant::BinSettings& ybins, const ant::BinSettings& zbins)
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
