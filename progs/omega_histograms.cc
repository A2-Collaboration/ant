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

#include <algorithm>

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


inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

void removesubstr(string& str, const string& sub) {

    string::size_type pos = 0;

    while(true) {

        pos = str.find(sub, pos);

        if(pos == str.npos)
            break;

        str.erase(pos,sub.length());
    }
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

    auto hists = GetListOf<TH1D>(dir, [] (const string& name) { return ends_with(name, "__gg");});

    hists.sort([] (const TH1* a, const TH1* b) { return a->GetEntries() > b->GetEntries();});

    hstack gg_stack("2#gamma IM");

    const std::vector<Color_t> cols = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
    auto cit = getCirculatIterator(cols.begin(), cols.end());

    int i=0;
    for(auto& h : hists) {

        h->SetFillColor(*(cit++));
        string title(h->GetTitle());
        removesubstr(title, "Pluto_dilepton ");
        h->SetTitle(title.c_str());

        gg_stack << h;
        if(i++ == 9)
            break;
    }

    canvas("gg im") << gg_stack << endc;

    app.Run();


    return 0;
}
