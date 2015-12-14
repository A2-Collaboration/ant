#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TRint.h"
#include "analysis/plot/root_draw.h"
#include <string>
#include <list>
#include "TH1D.h"
#include "base/CmdLine.h"
#include "TDirectory.h"
#include "TPad.h"
#include "TLine.h"
#include "base/interval.h"

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


class StringColorManager {
protected:
    std::map<std::string, Color_t> used;
    static const std::vector<Color_t> cols;

public:

    Color_t Get(const std::string& s) {
        const auto entry = used.find(s);
        if(entry == used.end()) {
            const auto index = used.size() % cols.size();
            const auto color = cols.at(index);
            used.emplace(s,color);
            return color;
        } else {
            return entry->second;
        }
    }

    void PrintList() const {
        for(const auto& entry : used) {
            cout << entry.first << endl;
        }
    }

};

const std::vector<Color_t> StringColorManager::cols = {kRed, kGreen+1, kBlue, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};

StringColorManager colors;

const static auto show_top_n = 5;

void DrawCutLines(TVirtualPad* pad, const ant::interval<double> x_cut) {
    TVirtualPad* p = gPad;
    pad->cd();
    double x1,y1,x2,y2;
    pad->GetRangeAxis(x1,y1,x2,y2);
    TLine* l1 = new TLine(x_cut.Start(),y1,x_cut.Start(),y2);
    l1->SetLineWidth(2);

    TLine* l2 = new TLine(x_cut.Stop(),y1,x_cut.Stop(),y2);
    l2->SetLineWidth(2);

    l1->Draw("same");
    l2->Draw("same");
    p->cd();
}

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

    colors.PrintList();

    app.Run();


    return 0;
}


void ggStack(TDirectory *dir) {
    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return std_ext::string_ends_with(name, "__gg");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack gg_stack("2#gamma IM", "2 #gamma invariant mass after #omega cut");

    int i=0;
        TH1D* sum = nullptr;
    for(auto& h : hists) {
        if(!sum) {
           sum = (TH1D*)h->Clone();
           sum->SetLineWidth(2);
           sum->SetLineColor(kBlack);
           sum->SetTitle("Sum");
           gg_stack << sum;
        } else {
            sum->Add(h);
        }
        string title(h->GetTitle());
        std_ext::removesubstr(title, "Pluto_dilepton ");
        std_ext::removesubstr(title, "2#gamma ");
        h->SetTitle(title.c_str());
        h->SetLineWidth(2);
        h->SetLineColor(colors.Get(title));

        gg_stack << h;
        if(i++ == show_top_n)
            break;
    }

    canvas c("gg im");
    c << padoption::Legend << drawoption("nostack")<< gg_stack << endc;
}

void mmStack(TDirectory *dir) {
    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return std_ext::string_ends_with(name, "__mm");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack gg_stack("MM");

    int i=0;
    TH1D* sum = nullptr;
    for(auto& h : hists) {

        if(!sum) {
           sum = (TH1D*)h->Clone();
           sum->SetLineWidth(2);
           sum->SetLineColor(kBlack);
           sum->SetTitle("Sum");
           gg_stack << sum;
        } else {
            sum->Add(h);
        }
        string title(h->GetTitle());
        std_ext::removesubstr(title, "Pluto_dilepton ");
        std_ext::removesubstr(title, "MM ");
        h->SetTitle(title.c_str());
        h->SetLineWidth(2);
        h->SetLineColor(colors.Get(title));

        gg_stack << h;
        if(i == show_top_n)
            break;
        ++i;
    }

    canvas c("mm");
    c << padoption::Legend << drawoption("nostack") << gg_stack << endc;

}

void gggStack(TDirectory *dir) {
    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return std_ext::string_ends_with(name, "__ggg");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack ggg_stack("3#gamma IM", "3 #gamma invariant mass");

    int i=0;
    TH1D* sum = nullptr;
    for(auto& h : hists) {

        if(!sum) {
           sum = (TH1D*)h->Clone();
           sum->SetLineWidth(2);
           sum->SetLineColor(kBlack);
           sum->SetTitle("Sum");
           ggg_stack << sum;
        } else {
            sum->Add(h);
        }
        string title(h->GetTitle());
        std_ext::removesubstr(title, "Pluto_dilepton ");
        std_ext::removesubstr(title, "3#gamma ");
        h->SetTitle(title.c_str());
        h->SetLineWidth(2);
        h->SetLineColor(colors.Get(title));

        ggg_stack << h;
        if(i++ == show_top_n)
            break;
    }

    canvas c("ggg im");
    c << padoption::Legend << drawoption("nostack") << ggg_stack << endc;
    DrawCutLines(gPad, interval<double>(680,780));
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
    c << padoption::Legend << drawoption("nostack") << ggg_stack << endc;
}
