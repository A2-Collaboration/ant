#include "HistogramFactory.h"

#include "base/std_ext/string.h"

#include "TDirectory.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace ant;
using namespace ant::analysis;
using namespace ant::std_ext;
using namespace std;


void HistogramFactory::goto_dir() const
{
    if(my_directory)
        my_directory->cd();
}

string HistogramFactory::MakeTitle(const string& title) const
{
    if(title_prefix.empty())
        return title;
    return std_ext::formatter() << title_prefix << ": " << title;
}

TDirectory *HistogramFactory::mkDirNumbered(const string &name, TDirectory *rootdir)
{
    TDirectory* dir = nullptr;

    unsigned n=0;
    do {
        const string dn = (n!=0) ? name+"_"+to_string(n) : name;
        ++n;

        rootdir->GetObject(dn.c_str(), dir);

        if(!dir) {

            dir = rootdir->mkdir(dn.c_str());

            if(!dir)
                throw("Can't create output directory \"" + dn +"\"");
        } else {
            dir = nullptr;
        }

    } while(dir==nullptr);

    return dir;
}

string HistogramFactory::GetNextName(const string &name, const string& autogenerate_prefix) const
{
    if(name.empty()) {
        if(!autogenerate_prefix.empty())
            return formatter() << autogenerate_prefix << setfill('0') << setw(3) << n_unnamed++;
        throw Exception("Cannot generate object with empty name in directory "+string(my_directory->GetPath()));
    } else {
        // check if name exists
        auto existing = my_directory->Get(name.c_str());
        if(existing != nullptr)
            throw Exception("Object with name "+name+" already exists in directory "+string(my_directory->GetPath()));
        return name;
    }
}

HistogramFactory::HistogramFactory(const string &directory_name, TDirectory* root, const string& title_prefix_):
    title_prefix(title_prefix_)
{

    if(!root)
        root=gDirectory;

    my_directory = mkDirNumbered(directory_name, root);

}


HistogramFactory::HistogramFactory(const string& directory_name, const HistogramFactory& parent, const string& title_prefix_)
    : HistogramFactory(directory_name, parent.my_directory,
                       // if the parent's title_prefix is empty, use the given one as is (possibly also empty)
                       parent.title_prefix.empty() ?
                           title_prefix_ :
                           // if the given prefix is empty, use the parent's title prefix (possibly also empty)
                           (title_prefix_.empty() ? parent.title_prefix :
                                                    // if both prefixes are not empty, construct one with colon in between
                                                    // this complicated procedure avoids inserting too many :, hopefully.
                                                    std_ext::formatter() << parent.title_prefix << ": " << title_prefix_))

{
}

void HistogramFactory::SetTitlePrefix(const string& title_prefix_)
{
    title_prefix = title_prefix_;
}

void HistogramFactory::SetDirDescription(const string &desc)
{
    my_directory->SetTitle(desc.c_str());
}

TH1D *HistogramFactory::makeTH1D(
        const string &title,
        const string &xlabel,
        const string &ylabel,
        const BinSettings& xbins,
        const string &name, bool sumw2) const
{
    return makeTH1D(title, {xlabel, xbins}, ylabel, name, sumw2);
}

TH1D*HistogramFactory::makeTH1D(
        const string& title,
        const AxisSettings& x_axis_settings,
        const string& name, bool sumw2) const
{
    // y axis denotes entries, so it's quite common that it's empty
    /// \todo think about automagically generating "per binwidth" label
    /// from given xbins in x_axis_settings
    return makeTH1D(title, x_axis_settings, "", name, sumw2);
}

TH1D* HistogramFactory::makeTH1D(
        const string& title,
        const AxisSettings& x_axis_settings,
        const string& ylabel, const string& name, bool sumw2) const
{
    auto& xbins = x_axis_settings;

    auto r = make<TH1D>(GetNextName(name).c_str(), MakeTitle(title).c_str(),
                        xbins.Bins(), xbins.Start(), xbins.Stop());

    r->SetXTitle(x_axis_settings.Label().c_str());
    r->SetYTitle(ylabel.c_str());

    if(sumw2) r->Sumw2();
    return r;
}

TH1D* HistogramFactory::makeTH1Dvarbin(
        const string& title,
        const vector<double>& edges,
        const string& xlabel, const string& ylabel,
        const string& name, bool sumw2) const
{
    auto r = make<TH1D>(GetNextName(name).c_str(), MakeTitle(title).c_str(),
                        edges.size()-1, &edges[0]);

    r->SetXTitle(xlabel.c_str());
    r->SetYTitle(ylabel.c_str());

    if(sumw2) r->Sumw2();
    return r;
}

TH2D* HistogramFactory::makeTH2D(
        const string &title,
        const string &xlabel,
        const string &ylabel,
        const BinSettings &xbins,
        const BinSettings &ybins,
        const string &name, bool  sumw2) const
{
    return makeTH2D(title, {xlabel, xbins}, {ylabel, ybins}, name, sumw2);
}

TH2D* HistogramFactory::makeTH2D(
        const string& title,
        const AxisSettings& x_axis_settings,
        const AxisSettings& y_axis_settings,
        const string& name, bool sumw2) const
{
    auto& xbins = x_axis_settings;
    auto& ybins = y_axis_settings;

    auto h = make<TH2D>(GetNextName(name).c_str(), MakeTitle(title).c_str(),
                         xbins.Bins(), xbins.Start(), xbins.Stop(),
                         ybins.Bins(), ybins.Start(), ybins.Stop());

    h->SetXTitle(x_axis_settings.Label().c_str());
    h->SetYTitle(y_axis_settings.Label().c_str());

    if(sumw2) h->Sumw2();
    return h;
}

TH3D* HistogramFactory::makeTH3D(
        const string &title,
        const string &xlabel,
        const string &ylabel,
        const string &zlabel,
        const BinSettings &xbins,
        const BinSettings &ybins,
        const BinSettings &zbins,
        const string& name, bool  sumw2) const
{
    return makeTH3D(title, {xlabel, xbins}, {ylabel, ybins}, {zlabel, zbins}, name, sumw2);
}

TH3D* HistogramFactory::makeTH3D(
        const string& title,
        const AxisSettings& x_axis_settings,
        const AxisSettings& y_axis_settings,
        const AxisSettings& z_axis_settings,
        const string& name, bool sumw2) const
{
    auto& xbins = x_axis_settings;
    auto& ybins = y_axis_settings;
    auto& zbins = z_axis_settings;

    auto h = make<TH3D>(GetNextName(name).c_str(), MakeTitle(title).c_str(),
                       xbins.Bins(), xbins.Start(), xbins.Stop(),
                       ybins.Bins(), ybins.Start(), ybins.Stop(),

                        zbins.Bins(), zbins.Start(), zbins.Stop());
    h->SetXTitle(x_axis_settings.Label().c_str());
    h->SetYTitle(y_axis_settings.Label().c_str());
    h->SetZTitle(z_axis_settings.Label().c_str());

    if(sumw2) h->Sumw2();
    return h;
}

TGraph* HistogramFactory::makeGraph(
        const string& title,
        const string& name) const
{
    auto g = new TGraph();

    g->SetName(GetNextName(name,"graph").c_str());
    g->SetTitle(title.c_str());
    g->SetMarkerStyle(kPlus);

    DirStackPush dirstack(*this);
    gDirectory->Add(g);
    return g;
}

TGraphErrors* HistogramFactory::makeGraphErrors(
        const string& title,
        const string& name) const
{
    auto g = new TGraphErrors();

    g->SetName(GetNextName(name,"graph").c_str());
    g->SetTitle(title.c_str());

    DirStackPush dirstack(*this);
    gDirectory->Add(g);

    return g;
}

TTree* HistogramFactory::makeTTree(const string& name) const
{
    // trees do not autogenerate histograms, thus provide empty prefix
    return make<TTree>(GetNextName(name, "").c_str(), MakeTitle(name.c_str()).c_str());
}

void HistogramFactory::addHistogram(TH1* h) const {
    h->SetDirectory(nullptr);
    DirStackPush dirstack(*this);
    gDirectory->Add(h);
}

HistogramFactory::DirStackPush::DirStackPush(const HistogramFactory& hf): dir(gDirectory)
{
    hf.goto_dir();
}

HistogramFactory::DirStackPush::~DirStackPush()
{
    dir->cd();
}
