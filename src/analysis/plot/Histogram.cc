#include "Histogram.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <sstream>
#include <iomanip>

using namespace ant::analysis;
using namespace std;


EColor HistogramFactory::GetNextColor()
{
    EColor c=kBlack;

    if(loopColors) {
        ++color;
        if( color == colors.end() )
            color = colors.begin();
        c=*color;
    }
    return c;
}

UInt_t HistogramFactory::GetNextHistnum()
{
    return histnum++;
}

HistogramFactory* HistogramFactory::instance = new HistogramFactory();


//TH1D* HistogramFactory::Make1D(const string &title, const string &xlabel, const string &ylabel, const BinSettings& bins, const string &name)


TH2D *HistogramFactory::Make2D(const std::string& title, const std::string& xlabel, const std::string& ylabel, const BinSettings& xbins, const BinSettings& ybins, const string &name)
{

    TH2D* h = new TH2D( GetNextHistName(name).c_str(), title.c_str(), xbins.Bins(), xbins.Start(), xbins.Stop(), ybins.Bins(), ybins.Start(), ybins.Stop());
    h->SetXTitle(xlabel.c_str());
    h->SetYTitle(ylabel.c_str());
    return h;
}

TH3D *HistogramFactory::Make3D(const string &title,
                               const string &xlabel,
                               const string &ylabel,
                               const string &zlabel,
                               const BinSettings &xbins,
                               const BinSettings &ybins,
                               const BinSettings &zbins,
                               const string &name)
{
    TH3D* h = new TH3D(
                GetNextHistName(name).c_str(),
                title.c_str(),
                xbins.Bins(), xbins.Start(), xbins.Stop(),
                ybins.Bins(), ybins.Start(), ybins.Stop(),
                zbins.Bins(), zbins.Start(), zbins.Stop());

    h->SetXTitle(xlabel.c_str());
    h->SetYTitle(ylabel.c_str());
    h->SetZTitle(zlabel.c_str());

    return h;
}

void HistogramFactory::ApplySettings(TH1D *hist, const string &title, const string &xlabel, const string &ylabel)
{
    if(!title.empty())
        hist->SetTitle(title.c_str());
    if(!xlabel.empty())
        hist->SetXTitle(xlabel.c_str());
    if(!ylabel.empty())
        hist->SetYTitle(ylabel.c_str());
    hist->SetLineColor(GetNextColor());
}

void HistogramFactory::ApplySettings(TH2D *hist, const string &title, const string &xlabel, const string &ylabel)
{
    if(!title.empty())
        hist->SetTitle(title.c_str());
    if(!xlabel.empty())
        hist->SetXTitle(xlabel.c_str());
    if(!ylabel.empty())
        hist->SetYTitle(ylabel.c_str());
}

string HistogramFactory::GetNextHistName(const std::string& name)
{
    stringstream s;

    if(name.empty()) {
        s << "hist" << setfill('0') << setw(3) << GetNextHistnum();
    } else {
        s << name;
    }
    return s.str();
}

HistogramFactory::HistogramFactory():
    color(colors.begin()),
    histnum(0),
    loopColors(false)
{
}

void HistogramFactory::ResetColors()
{
    color=colors.begin();
}

// Color set for 1D histogram lines.
std::vector<EColor> HistogramFactory::colors = {kBlue, kRed, kGreen, kMagenta, kCyan, kBlack};
