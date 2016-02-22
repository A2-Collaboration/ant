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

#include "base/std_ext/iterators.h"
#include "base/std_ext/string.h"

#include "root-addons/analysis_codes/hstack.h"

using namespace ant;
using namespace std;


template<class T, typename Func>
void GetListOf(TDirectory* dir,
               std::list<T*>& list,
               Func match=[] (const string&, size_t){return true;},
               size_t level = 0)
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
            GetListOf(dynamic_cast<TDirectoryFile*>(key->ReadObj()), list, match, level+1);
            continue;
        }

        const auto& name = key->GetName();
        if(match(name, level)) {
            objectPtr = dynamic_cast<T*>(key->ReadObj());
            if ( !objectPtr )
                continue;
            list.push_back(objectPtr);
        }
    }
}

void do_Sum_Signal_Bkg(TDirectory* dir,
                       const string& hname, const string& htitle,
                       const unsigned rebin = 1,
                       const bool bkgdetail = false);

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

    void Clear() {
        used.clear();
    }

    void PrintList() const {
        for(const auto& entry : used) {
            cout << entry.first << endl;
        }
    }

};

const std::vector<Color_t> StringColorManager::cols = {kBlack, kRed, kGreen+1, kBlue, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray+1, kTeal, kGray};

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

    TDirectory* dir = nullptr;

    dir = file->GetDirectory("EtapOmegaG/Sig");
    if(dir == nullptr) {
        LOG(ERROR) << "Directory not found in file";
        return 1;
    }


    do_Sum_Signal_Bkg(dir, "Chi2_Min", "#chi^{2} Minimum", 4);
    do_Sum_Signal_Bkg(dir, "g_EtaPrime_E", "E_{#gamma} in #eta' frame", 4);
    do_Sum_Signal_Bkg(dir, "g_EtaPrime_E", "E_{#gamma} in #eta' frame", 4, true);
    do_Sum_Signal_Bkg(dir, "steps", "Cut efficiencies", 1, true);

    colors.Clear();
    dir = file->GetDirectory("EtapOmegaG/Ref");
    if(dir == nullptr) {
        LOG(ERROR) << "Directory not found in file";
        return 1;
    }
    do_Sum_Signal_Bkg(dir, "IM_etap", "2#gamma IM", 2);
    do_Sum_Signal_Bkg(dir, "IM_etap", "2#gamma IM", 2, true);

    colors.PrintList();

    app.Run();

    return 0;
}

TH1D* MakeSum(std::list<TH1D*>::const_iterator start,
              std::list<TH1D*>::const_iterator stop,
              const std::string& title = "Sum")
{
    TH1D* sum = reinterpret_cast<TH1D*>((*(start++))->Clone());
    sum->SetTitle(title.c_str());
    sum->SetLineColor(colors.Get(sum->GetTitle()));
    for(auto it = start; it != stop; it++)
        sum->Add(*it);
    return sum;
}

void sanitize_title(TH1D* h) {
    string title(h->GetTitle());
    title = title.substr(0, title.find(':'));
    h->SetTitle(title.c_str());
}

void do_Sum_Signal_Bkg(TDirectory *dir,
                       const std::string& hname,
                       const std::string& htitle,
                       const unsigned rebin,
                       const bool bkgdetail
                       ) {
    std::list<TH1D*> hists;
    GetListOf(dir, hists, [hname] (const string& name, size_t level) { return name == hname && level>0;});
    for(TH1D* h : hists) {
        if(rebin>1)
            h->Rebin(rebin);
        h->SetLineWidth(2);
    }

    hstack mystack(hname,htitle);

    // assume that hists are ordered, first one is signal,
    // the rest are background channels
    mystack << MakeSum(hists.begin(), hists.end(), "Sum");
    mystack << MakeSum(hists.begin(), std::next(hists.begin()), "Signal");
    if(bkgdetail) {
        for(auto it = std::next(hists.begin()); it != hists.end(); it++) {
            TH1D* h = *it;
            sanitize_title(h);
            h->SetLineColor(colors.Get(h->GetTitle()));
            mystack << h;
        }
    }
    else {
        mystack << MakeSum(std::next(hists.begin()), hists.end(), "Background");
    }

    canvas c("c_"+hname+(bkgdetail ? "_detail" : ""));
    c << padoption::Legend << drawoption("nostack") << &mystack << endc;
}



