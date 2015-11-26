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
void GetListOf(TDirectory* dir,
               std::list<T*>& list,
               Func match=[] (const string&){return true;})
{
    TList* keys = dir->GetListOfKeys();

    if(!keys)
        return;

    T* objectPtr = nullptr;
    TKey* key = nullptr;
    TIter nextk(keys);

    while((key = (TKey*)nextk()))
    {
        const std::string classname(key->GetClassName());

        if(classname == "TDirectoryFile") {
            GetListOf(dynamic_cast<TDirectoryFile*>(key->ReadObj()), list, match);
            continue;
        }

        const auto& name = key->GetName();
        if(match(name)) {
            objectPtr = dynamic_cast<T*>(key->ReadObj());
            if ( !objectPtr )
                continue;
            list.push_back(objectPtr);
        }
    }
}

void do_stack(TDirectory* dir);

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

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");

    auto cmd_file = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");;

    cmd.parse(argc, argv);

    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("plot",&fake_argc,fake_argv,nullptr,0,true);

    TFile* file = new TFile(cmd_file->getValue().c_str());

    TDirectory* dir = file->GetDirectory("EtapOmegaG/Sig");
    if(dir == nullptr) {
        LOG(ERROR) << "Directory not found in file";
        return 1;
    }


    do_stack(dir);

    colors.PrintList();

    app.Run();

    return 0;
}


void do_stack(TDirectory *dir) {
    std::list<TH1D*> hists;
    GetListOf(dir, hists, [] (const string& name) { return name == "Chi2_Min";});

    //hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    //hists.resize(5);

    hstack gg_stack("2#gamma IM", "2 #gamma invariant mass after #omega cut");

    TH1D* sum = nullptr;
    for(auto& h : hists) {
        h->Rebin(4);
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
    }

    canvas c("gg im");
    c << padoption::set(padoption_t::Legend) << drawoption("nostack")<< gg_stack << endc;
}



