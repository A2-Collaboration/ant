#include "TreeTools.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCut.h"
#include "TDirectory.h"

using namespace std;

TH1*Draw(TTree* tree, const string& formula, const TCut& cut, const int bins, const double min, const double max) {
    static unsigned n = 0;
    const char* hname = Form("h_%d", n++);

    tree->Draw(Form("%s>>%s(%d,%lf,%lf)",formula.c_str(),hname,bins,min,max),cut);
    TH1* h = NULL;
    gDirectory->GetObject(hname, h);
    return h;
}

TH2*Draw(TTree* tree, const string& formula, const TCut& cut, const int xbins, const double xmin, const double xmax, const int ybins, const double ymin, const double ymax) {
    static unsigned n = 0;
    const char* hname = Form("h_%d", n++);

    tree->Draw(Form("%s>>%s(%d,%lf,%lf,%d,%lf,%lf)",formula.c_str(),hname,xbins,xmin,xmax,ybins,ymin,ymax),cut,"colz");
    TH2* h = NULL;
    gDirectory->GetObject(hname, h);
    return h;
}
