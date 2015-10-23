#include "plot/HistogramFactories.h"

#include "plot/SmartHist.h"
#include "data/Particle.h"

#include "TMath.h"
#include "TDirectory.h"
#include "TTree.h"

using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace std;


TDirectory *SmartHistFactory::goto_dir()
{
    TDirectory* old = gDirectory;
    if(dir)
        dir->cd();

    return old;
}

void SmartHistFactory::restore_dir(TDirectory* dir)
{
    if(dir)
        dir->cd();
}

string SmartHistFactory::MakeTitle(const string& title)
{
    if(title_prefix.empty())
        return title;
    return title_prefix+": "+title;
}

SmartHistFactory::SmartHistFactory(const string &directory_name, TDirectory* root, const string& title_prefix_)
  : dir(), base_factory(), title_prefix(title_prefix_)
{

    if(!root)
        root=gDirectory;

    dir = root->mkdir(directory_name.c_str());

    if(!dir)
        dir=gDirectory;
}


SmartHistFactory::SmartHistFactory(const string &directory_name, const SmartHistFactory& parent, const string& title_prefix_)
  : dir(), base_factory(), title_prefix(title_prefix_)
{
    dir = parent.dir->mkdir(directory_name.c_str());
    if(!dir)
        dir=gDirectory;
}

void SmartHistFactory::SetTitlePrefix(const string& title_prefix_)
{
    title_prefix = title_prefix_;
}




TH2D *SmartHistFactory::makeTH2D(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings &ybins, const string &name)
{
    TDirectory* old = goto_dir();
    TH2D* r = base_factory.Make2D(MakeTitle(title),xlabel,ylabel,xbins,ybins,name);
    restore_dir(old);
    return r;
}

TH3D *SmartHistFactory::makeTH3D(const string &title,
                                 const string &xlabel,
                                 const string &ylabel,
                                 const string &zlabel,
                                 const BinSettings &xbins,
                                 const BinSettings &ybins,
                                 const BinSettings &zbins,
                                 const string &name)
{
    TDirectory* old = goto_dir();
    TH3D* r = base_factory.Make3D(MakeTitle(title),xlabel,ylabel,zlabel,xbins,ybins,zbins,name);
    restore_dir(old);
    return r;
}

TTree*SmartHistFactory::makeTTree(const string& name)
{
    TDirectory* old = goto_dir();
    TTree* t = new TTree(name.c_str(), name.c_str());
    restore_dir(old);
    return t;
}

SmartHist1<const ParticlePtr &> SmartHistFactory::InvariantMass(const string &title, const string &xlabel, const string &ylabel, BinSettings bins, const string &name)
{
    return makeHist<const ParticlePtr&>(
                [] (const ParticlePtr& p) { return p->M();},
            MakeTitle(title),
            xlabel,
            ylabel,
            bins,
    name);
}

SmartHist1<const ParticlePtr &> SmartHistFactory::ThetaAngle(const string &title, const string &xlabel, const string &ylabel, BinSettings bins, const string &name)
{
    return makeHist<const ParticlePtr&>(
                [] (const ParticlePtr& p) { return p->Theta() * TMath::DegToRad();},
            MakeTitle(title),
            xlabel,
            ylabel,
            bins,
            name);
}

SmartHist1<const ParticlePtr &> SmartHistFactory::KinEnergyPlot(const string &title, const string &xlabel, const string &ylabel, BinSettings bins, const string &name)
{
    return makeHist<const ParticlePtr&>(
                [] (const ParticlePtr& p) { return p->Ek();},
            MakeTitle(title),
            xlabel,
            ylabel,
            bins,
            name);
}

TH1D *SmartHistFactory::copyTH1D(TH1D *hist, const std::string& newname)
{
    TDirectory* old = goto_dir();
    TH1D* newhist = (TH1D*)hist->Clone(newname.c_str());
    restore_dir(old);
    return newhist;
}

//void SmartHistFactory::Add(SmartHist1Base &hist, const std::string& as_name)
//{
//    TDirectory* old = begin_make_histogram();
//        hist.histogram->SetName(as_name.c_str());
//        dir->Add(hist.histogram, true);

//        hist.cleanup = false;
//    end_make_histogram(old);
//}

//void SmartHistFactory::Add(SmartHist1Base &hist)
//{
//    Add(hist, base_factory.GetNextHistName());
//}
