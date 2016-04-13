#include "plot/HistogramFactories.h"

#include "base/std_ext/string.h"

#include "TMath.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <iomanip>

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
        const string dn = (n!=0) ? name+to_string(n) : name;
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

string HistogramFactory::GetNextHistName(const string &name) const
{
    if(name.empty()) {
        return formatter() << "hist" << setfill('0') << setw(3) << n_unnamed++;
    } else
        return name;
}

HistogramFactory::HistogramFactory(const string &directory_name, TDirectory* root, const string& title_prefix_):
    title_prefix(title_prefix_)
{

    if(!root)
        root=gDirectory;

    my_directory = mkDirNumbered(directory_name, root);

}


HistogramFactory::HistogramFactory(const string& directory_name, const HistogramFactory& parent, const string& title_prefix_)
  : my_directory(),
    title_prefix
    (
        parent.title_prefix.empty() ?
            title_prefix_ :
            (title_prefix_.empty() ? parent.title_prefix :
                                     std_ext::formatter() << parent.title_prefix << ": " << title_prefix_)
    )
{
    my_directory = parent.my_directory->mkdir(directory_name.c_str());
    if(!my_directory)
        my_directory=gDirectory;
}

void HistogramFactory::SetTitlePrefix(const string& title_prefix_)
{
    title_prefix = title_prefix_;
}

string HistogramFactory::GetTitlePrefix() const
{
    return title_prefix;
}

void HistogramFactory::SetDirDescription(const string &desc)
{
    my_directory->SetTitle(desc.c_str());
}

TH1D *HistogramFactory::makeTH1D(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name) const
{
    auto r = make<TH1D>(GetNextHistName(name).c_str(), MakeTitle(title).c_str(),
                        bins.Bins(), bins.Start(), bins.Stop());
    r->SetXTitle(xlabel.c_str());
    r->SetYTitle(ylabel.c_str());
    return r;
}




TH2D *HistogramFactory::makeTH2D(const string &title,
                                 const string &xlabel,
                                 const string &ylabel,
                                 const BinSettings &xbins,
                                 const BinSettings &ybins,
                                 const string &name) const
{
    auto h = make<TH2D>(GetNextHistName(name).c_str(), MakeTitle(title).c_str(),
                         xbins.Bins(), xbins.Start(), xbins.Stop(),
                         ybins.Bins(), ybins.Start(), ybins.Stop());
    h->SetXTitle(xlabel.c_str());
    h->SetYTitle(ylabel.c_str());
    return h;
}

TH3D *HistogramFactory::makeTH3D(const string &title,
                                 const string &xlabel,
                                 const string &ylabel,
                                 const string &zlabel,
                                 const BinSettings &xbins,
                                 const BinSettings &ybins,
                                 const BinSettings &zbins,
                                 const string &name) const
{
    auto h = make<TH3D>(GetNextHistName(name).c_str(), MakeTitle(title).c_str(),
                       xbins.Bins(), xbins.Start(), xbins.Stop(),
                       ybins.Bins(), ybins.Start(), ybins.Stop(),
                       zbins.Bins(), zbins.Start(), zbins.Stop());
    h->SetXTitle(xlabel.c_str());
    h->SetYTitle(ylabel.c_str());
    h->SetZTitle(zlabel.c_str());
    return h;
}

TTree*HistogramFactory::makeTTree(const string& name)
{
    TTree* t = make<TTree>(name.c_str(), MakeTitle(name.c_str()).c_str());
    return t;
}

HistogramFactory::DirStackPush::DirStackPush(): dir(gDirectory)
{}

HistogramFactory::DirStackPush::~DirStackPush()
{
    dir->cd();
}
