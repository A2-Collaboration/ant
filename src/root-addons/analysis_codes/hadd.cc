#include "hadd.h"

#include "hstack.h"
#include "tree/TAntHeader.h"
#include "base/ProgressCounter.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TFileMergeInfo.h"

#include <algorithm>

using namespace std;
using namespace ant;


template<typename T>
struct pair_t {
    explicit pair_t(const string& name) : Name(name) {}
    string Name;
    T Item;
};

template<typename C, typename... Args>
void add_by_name(C& c, const string& name, Args&&... args) {
    using T = typename C::value_type;
    auto it = std::find_if(c.begin(), c.end(), [name] (const T& item) {
        return item.Name == name;
    });
    if(it == c.end()) {
        c.emplace_back(name);
        c.back().Item.emplace_back(std::forward<Args>(args)...);
    }
    else {
        it->Item.emplace_back(std::forward<Args>(args)...);
    }
}

void hadd::MergeRecursive(TDirectory& target, const hadd::sources_t& sources, unsigned& nPaths)
{

        nPaths++;
        ProgressCounter::Tick();

        vector<pair_t<sources_t>> dirs;

        vector<pair_t<unique_ptrs_t<TH1>>>    hists;
        vector<pair_t<unique_ptrs_t<hstack>>> stacks;
        vector<pair_t<unique_ptrs_t<TAntHeader>>> headers;

        for(auto& source : sources) {
            TList* keys = source->GetListOfKeys();
            if(!keys)
                continue;

            // first create a unique list of names,
            // this prevents object with different cycles
            TIter nextk(keys);
            string prev_keyname;
            while(auto key = dynamic_cast<TKey*>(nextk()))
            {
                const string keyname = key->GetName();
                if(prev_keyname == keyname)
                    continue;
                prev_keyname = keyname;

                auto cl = TClass::GetClass(key->GetClassName());

                if(cl->InheritsFrom(TDirectory::Class())) {
                    auto dir = dynamic_cast<TDirectory*>(key->ReadObj());
                    add_by_name(dirs, keyname, dir);
                }
                else if(cl->InheritsFrom(TH1::Class())) {
                    auto obj = dynamic_cast<TH1*>(key->ReadObj());
                    add_by_name(hists, keyname, obj);
                }
                else if(cl->InheritsFrom(hstack::Class())) {
                    auto obj = dynamic_cast<hstack*>(key->ReadObj());
                    add_by_name(stacks, keyname, obj);
                }
                else if(cl->InheritsFrom(TAntHeader::Class())) {
                    auto obj = dynamic_cast<TAntHeader*>(key->ReadObj());
                    add_by_name(headers, keyname, obj);
                }
            }
        }



        for(const auto& it_dirs : dirs) {
            auto newdir = target.mkdir(it_dirs.Name.c_str());
            MergeRecursive(*newdir, it_dirs.Item, nPaths);
        }

        target.cd();
        TFileMergeInfo info(addressof(target)); // for calling Merge

        for(const auto& it_hists : hists) {
            auto& items = it_hists.Item;

            // check if at least one hist has labels,
            // the others could be never filled (so ROOT treats them as normal hists)
            const auto hasLabels = [] (const unique_ptr<TH1>& h) {
                return h->GetXaxis()->GetLabels() != nullptr;
            };
            const auto it_h_withLabels = std::find_if(items.begin(), items.end(), hasLabels);

            auto& first = items.front();
            if(it_h_withLabels != items.end()) {

                // again, scan the histograms for empty hists without labels
                // IMHO, this is a bug in ROOT that empty hists cannot be merged with labeled hists
                auto& h_withLabels = *it_h_withLabels;
                for(auto& h : items) {
                    if(hasLabels(h))
                        continue;
                    for(int bin=0;bin<h->GetNbinsX()+1;bin++)
                        if(h->GetBinContent(bin) != 0)
                            throw std::runtime_error("Found non-empty unlabeled hist "
                                                     + string(h->GetDirectory()->GetPath()));
                    // prepare the axis labels of the empty hist, labeled hist should have at least
                    // one bin filled
                    h->Fill(h_withLabels->GetXaxis()->GetBinLabel(1), 0.0);
                }

                auto& first = items.front();
                TList c;
                for(auto it = next(items.begin()); it != items.end(); ++it) {
                    c.Add(it->get());
                }
                first->Merge(addressof(c));
            }
            else {
                for(auto it = next(items.begin()); it != items.end(); ++it) {
                    first->Add(it->get());
                }
            }
            target.WriteTObject(first.get());
        }

        for(const auto& it : stacks) {
            auto& items = it.Item;
            auto& first = items.front();
            TList c;
            for(auto it = next(items.begin()); it != items.end(); ++it) {
                c.Add(it->get());
            }
            first->Merge(addressof(c), addressof(info));
            target.WriteTObject(first.get());
        }

        for(const auto& it : headers) {
            auto& items = it.Item;
            auto& first = items.front();
            TList c;
            for(auto it = next(items.begin()); it != items.end(); ++it) {
                c.Add(it->get());
            }
            first->Merge(addressof(c));
            target.WriteTObject(first.get());
        }

}
