#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TRint.h"
#include "analysis/plot/root_draw.h"
#include <string>
#include <list>
#include "TH1D.h"
#include "base/CmdLine.h"
#include "TDirectory.h"

#include "base/iterators.h"
#include "base/std_ext/string.h"

using namespace ant;
using namespace std;


template<class T, typename Func>
std::list<T*> GetListOf(TDirectory* dir, Func func=[] (const string& name){return true;})
{
    std::list<T*> theList;

    TList* keys = dir->GetListOfKeys();

    if(!keys)
        return {};

    T* objectPtr = nullptr;
    TKey* key = nullptr;
    TIter nextk(keys);

    while((key = (TKey*)nextk()))
    {
        const auto& name = key->GetName();
        if(func(name)) {
            objectPtr = dynamic_cast<T*>(key->ReadObj());
            if ( !objectPtr )
                continue;
            theList.push_back(objectPtr);
        }
    }
    return theList;
}

void ggStack(TDirectory* dir);
void mmStack(TDirectory* dir);
void gggStack(TDirectory *dir);
void nPhotons(TDirectory* dir);

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("Ant", ' ', "0.1");

    auto cmd_files = cmd.add<TCLAP::MultiArg<string>>("i","input","Input files",true,"filename");

    cmd.parse(argc, argv);


    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("Ant",&fake_argc,fake_argv,nullptr,0,true);


    WrapTFileInput infiles;

    for(const auto& f : cmd_files->getValue()) {
        infiles.OpenFile(f);
    }

    TDirectory* dir = nullptr;

    infiles.GetObject("OmegaEtaG",dir);


    ggStack(dir);
    mmStack(dir);
    gggStack(dir);
    nPhotons(dir);

    app.Run();


    return 0;
}


void ggStack(TDirectory *dir) {
    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return std_ext::string_ends_with(name, "__gg");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack gg_stack("2#gamma IM");

    const std::vector<Color_t> cols = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
    auto cit = getCirculatIterator(cols.begin(), cols.end());

    int i=0;
    for(auto& h : hists) {

        h->SetFillColor(*(cit++));
        string title(h->GetTitle());
        std_ext::removesubstr(title, "Pluto_dilepton ");
        std_ext::removesubstr(title, "2#gamma ");
        h->SetTitle(title.c_str());

        gg_stack << h;
        if(i++ == 9)
            break;
    }

    canvas c("gg im");
    c << padoption::set(padoption_t::Legend) << gg_stack << endc;
}

void mmStack(TDirectory *dir) {
    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return std_ext::string_ends_with(name, "__mm");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack gg_stack("MM");

    const std::vector<Color_t> cols = {kRed, kGreen, kBlue, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
    auto cit = getCirculatIterator(cols.begin(), cols.end());

    int i=0;
    for(auto& h : hists) {

        h->SetLineColor(*(cit++));
        h->SetLineWidth(2);
        string title(h->GetTitle());
        std_ext::removesubstr(title, "Pluto_dilepton ");
        std_ext::removesubstr(title, "MM ");
        h->SetTitle(title.c_str());

        gg_stack << h;
        if(i == 3)
            break;
        ++i;
    }

    canvas c("mm");
    c << padoption::set(padoption_t::Legend) << drawoption("nostack") << gg_stack << endc;
}

void gggStack(TDirectory *dir) {
    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return std_ext::string_ends_with(name, "__ggg");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack ggg_stack("3#gamma IM");

    const std::vector<Color_t> cols = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
    auto cit = getCirculatIterator(cols.begin(), cols.end());

    int i=0;
    for(auto& h : hists) {

        h->SetFillColor(*(cit++));
        string title(h->GetTitle());
        std_ext::removesubstr(title, "Pluto_dilepton ");
        std_ext::removesubstr(title, "3#gamma ");
        h->SetTitle(title.c_str());

        ggg_stack << h;
        if(i++ == 9)
            break;
    }

    canvas c("ggg im");
    c << padoption::set(padoption_t::Legend) << ggg_stack << endc;
}

void nPhotons(TDirectory *dir) {
    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return std_ext::string_ends_with(name, "__ncand");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack ggg_stack("ncand");

    const std::vector<Color_t> cols = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
    auto cit = getCirculatIterator(cols.begin(), cols.end());

    int i=0;
    for(auto& h : hists) {

        h->SetLineColor(*(cit++));
        string title(h->GetTitle());
        std_ext::removesubstr(title, "Pluto_dilepton ");
        std_ext::removesubstr(title, ": # Photons");
        h->SetTitle(title.c_str());

        ggg_stack << h;
        if(i++ == 9)
            break;
    }

    canvas c("n Photons");
    c << padoption::set(padoption_t::Legend) << drawoption("nostack") << ggg_stack << endc;
}
