#include "hstack.h"



#include "base/std_ext/string.h"
#include "base/std_ext/math.h"
#include "base/Logger.h"


#include "base/cereal/cereal.hpp"
#include "base/cereal/types/string.hpp"
#include "base/cereal/types/vector.hpp"
#include "base/cereal/archives/binary.hpp"
#include "tree/stream_TBuffer.h"

#include "TBrowser.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TH1.h"
#include "THStack.h"
#include "TDirectory.h"

#include <list>

using namespace ant;
using namespace std;

const hstack::options_t hstack::options_t::all_enabled = {true, true, true, true};

double hstack::Global_MC_Scaling = std_ext::NaN;
map<TH1*, double> hstack::Scaled_Hists = {};

interval<double> hstack::GlobalYAxisRange = {std_ext::NaN, std_ext::NaN};
interval<interval<double>> hstack::GlobalLegendPosition = {
    {0.7,  0.63},
    {0.88, 0.88}
};


bool hstack::options_t::operator==(const hstack::options_t& rhs) const {
    /// \todo check if one could use the serialize method here...
    auto make_tuple = [] (const options_t& o) {
        return std::tie(o.UseIntelliLegend, o.IgnoreEmptyHist,
                        o.DrawNoStack, o.ShowEntriesInLegend);
    };
    return make_tuple(*this) == make_tuple(rhs);
}

hstack::hstack(const string& name, const std::string& title,
               const options_t& options_) :
    TNamed(name.c_str(), title.c_str()),
    options(options_)
{
    gDirectory->Append(this);
}

hstack::hstack() {}

hstack::~hstack()
{}

TH1* hstack::hist_t::GetPtr(const string& path)
{
    auto ptr = dynamic_cast<TH1*>(gDirectory->Get(path.c_str()));
    if(ptr == nullptr)
        LOG(WARNING) << "Could not Get() TH1* for path " << path;
    return ptr;
}

string hstack::hist_t::GetPath(const TH1* ptr)
{
    const std::string path = std_ext::formatter()
                             << ptr->GetDirectory()->GetPath()
                             << "/" << ptr->GetName();
    // strip the file of the full path
    const auto n = path.find_first_of(':');
    if(n == path.npos)
        return path;
    return path.substr(n+1);
}

hstack& hstack::operator<<(TH1* hist)
{

    hists.emplace_back(hist, current_option);

    xlabel = hist->GetXaxis()->GetTitle();
    ylabel = hist->GetYaxis()->GetTitle();

    current_option.Z = 0; // reset z after each add

    return *this;
}

hstack& hstack::operator<<(const drawoption& c)
{
    current_option.DrawOption = c.Option();
    return *this;
}

hstack& hstack::operator<<(const ModOption_t& option)
{
    current_option = option;
    return *this;
}

ostream& hstack::Print(ostream& s) const
{
    s << "ant::hstack: "  << GetName() << ": " << GetTitle() << "\n";
    if(isfinite(Global_MC_Scaling))
        s << "Global MC Scaling = " << Global_MC_Scaling;
    /// \todo print current options...?
    for(const auto& h : hists) {
        s << "  " << h.Path << "\n";
    }
    s << endl;
    return s;
}



void hstack::Print(const char*) const
{
    cout << *this;
}

void hstack::Print() const
{
    Print("");
}

void hstack::checkHists()
{
    auto it_hist = hists.begin();
    while(it_hist != hists.end()) {
        hist_t hist = *it_hist;
        if(hist.Ptr == nullptr) {
            // try to get it now, might be that after merge
            // the Ptr is not yet initialized
            hist.Ptr = hist_t::GetPtr(hist.Path);
            // erase if still unsuccessful
            if(hist.Ptr == nullptr) {
                it_hist = hists.erase(it_hist);
                continue;
            }
        }
        it_hist++;
    }
}



bool hstack::IsCompatible(const hstack& other) const
{
    // we do not check hists here, since that's what could be merged
    return string(GetName()) == string(other.GetName()) &&
            string(GetTitle()) == string(other.GetTitle()) &&
            xlabel == other.xlabel &&
            ylabel == other.ylabel &&
            options == other.options;

}

template<typename Att>
void AttCopy(const Att* src, Att* dest) {
    src->Copy(*dest);
}

void hstack::buildIntelliLegend() const
{
    const auto& p = GlobalLegendPosition;
    auto legend = new TLegend(p.Start().Start(), // x1
                              p.Start().Stop(),  // y1
                              p.Stop().Start(),  // x2
                              p.Stop().Stop()    // y2
                              );

    const auto& delim = ": "; // this makes the tokens

    vector<vector<string>> title_parts;
    map<string, unsigned>  token_counter;
    for(const auto& hist : hists) {
        const string title = hist.Ptr->GetTitle();
        title_parts.emplace_back(std_ext::tokenize_string(title, delim));
        for(const auto& token : title_parts.back())
            token_counter[token]++;
    }

    vector<vector<string>> unique_title_parts;
    for(size_t i=0; i< title_parts.size(); i++) {
        const hist_t& hist = hists[i];
        if(options.IgnoreEmptyHist && hist.Ptr->GetEntries()==0)
            continue;

        const auto& tokens = title_parts[i];
        unique_title_parts.emplace_back();
        auto& unique_tokens = unique_title_parts.back();
        for(const auto& token : tokens)
            if(token_counter[token] < hists.size())
                unique_tokens.emplace_back(token);
        string unique_title = std_ext::concatenate_string(unique_tokens, delim);
        if(options.ShowEntriesInLegend)
            unique_title += std_ext::formatter() << " (" << hists[i].Ptr->GetEntries() << ")";

        auto entry = legend->AddEntry((TObject*)0, unique_title.c_str());
        AttCopy<TAttLine>(hist.Ptr, entry);
        AttCopy<TAttFill>(hist.Ptr, entry);
        AttCopy<TAttMarker>(hist.Ptr, entry);
        if(hist.Option.BkgColor>=0) {
            // during legend building, modify fill color
            entry->SetFillColor(hist.Option.BkgColor);
            entry->SetFillStyle(1001);
            entry->SetLineWidth(2);
        }
    }
    legend->Draw();
}

void hstack::Draw(const char* option)
{
    if(hists.empty())
        return;

    // ensure the histograms are ok
    checkHists();

    // hack the THStack such that the last histogram added gets drawn first
    struct THStack_hack : THStack {
        using THStack::THStack; // use constructors
        virtual void Paint(const char* chopt="") override {

            if(GlobalYAxisRange.IsSane()) {
                SetMinimum(GlobalYAxisRange.Start());
                SetMaximum(GlobalYAxisRange.Stop());
            }
            else {
                // resets it to default -1111, well, THStack knows what to do then
                SetMinimum();
                SetMaximum();
            }

            if(fHists->IsEmpty())
                return;



            // this helper extracts the z number from the options,
            // and temporarily removes them from the the option before
            // painting the histogram
            struct tmphist_t {
                tmphist_t(TObject* obj, const char* option) :
                    Obj(obj),
                    OrigOption(option)
                {
                    stringstream ss_option;
                    ss_option << OrigOption;
                    // try to extract Z (see hstack->Add as well)
                    ss_option >> Z;
                    // use rest as is
                    ss_option >> Option;
                }
                TObject* Obj;
                string OrigOption;
                string Option;
                int Z = 0;
                bool operator< (const tmphist_t& other) const {
                    return Z < other.Z;
                }
            };

            list<tmphist_t> tmp_hists;
            // we need to iterate over links, since we need the option string
            // attached to each item
            auto firstlink = fHists->FirstLink();
            while(firstlink) {
                tmp_hists.emplace_back(firstlink->GetObject(), firstlink->GetOption());
                firstlink = firstlink->Next();
            }
            // save for later
            auto tmp_hists_backup = tmp_hists;

            // reverse first
            std::reverse(tmp_hists.begin(), tmp_hists.end());

            // then sort according to Z (the sort is stable)
            tmp_hists.sort();

            fHists->Clear("nodelete");
            for(const auto& o : tmp_hists) {
                fHists->Add(o.Obj, o.Option.c_str());
            }

            // let's hope that this does not throw an exception
            // before we unreverse the fHists
            THStack::Paint(chopt);

            // restore from backup
            fHists->Clear("nodelete");
            for(const auto& o : tmp_hists_backup) {
                fHists->Add(o.Obj, o.OrigOption.c_str());
            }

        }

    };

    auto stack = new THStack_hack((string(GetName())+"_").c_str(), GetTitle());

    unsigned nAdded = 0;
    for(const auto& hist : hists) {
        if(options.IgnoreEmptyHist && hist.Ptr->GetEntries()==0)
            continue;
        stringstream ss_option;
        ss_option << hist.Option.Z << hist.Option.DrawOption;
        stack->Add(hist.Ptr, ss_option.str().c_str());
        nAdded++;
    }

    if(nAdded==0) {
        LOG(WARNING) << "No histograms in ant::hstack to draw (maybe all empty)";
        return;
    }

    UpdateMCScaling();

    // finally draw the stack
    string option_str(option);
    if(options.DrawNoStack)
        option_str += "nostack";
    stack->Draw(option_str.c_str());

    // axis business
    auto xaxis = stack->GetXaxis();
    xaxis->SetTitle(xlabel.c_str());

    auto yaxis = stack->GetYaxis();
    yaxis->SetTitle(ylabel.c_str());

    if(options.UseIntelliLegend)
        buildIntelliLegend();
}

void hstack::Browse(TBrowser* b)
{
    Draw(b ? b->GetDrawOption() : "");
    gPad->Update();
}

void hstack::SetGlobalMCScaling(double scaling)
{
    if(isfinite(scaling) && scaling != 0) {
        Global_MC_Scaling = scaling;
        UpdateMCScaling();
        gPad->Modified();
        gPad->Update();
        cout << endl << "ant::hstack: Set global MC scaling to " << scaling << endl;
    }
}

void hstack::SetGlobalYAxisRange(double low, double high)
{
    GlobalYAxisRange = {low, high};
    // even if axis range is not sane, redraw
    gPad->Modified();
    gPad->Update();
}

void hstack::SetGlobalLegendPosition(double x1, double y1, double x2, double y2)
{
    GlobalLegendPosition = {{x1, y1}, {x2, y2}};
}

void hstack::UpdateMCScaling()
{
    if(isfinite(Global_MC_Scaling)) {
        for(hist_t& hist : hists) {
            /// \todo find better way to detect which histograms
            /// should NOT be scaled. We assume that data is always drawn
            /// with error bars.
            if(std_ext::contains(hist.Option.DrawOption, "E"))
                continue;
            // have a look if this hist was scaled already
            double scale = Global_MC_Scaling;
            auto it_scaled_hist = Scaled_Hists.find(hist.Ptr);
            if(it_scaled_hist != Scaled_Hists.end()) {
                scale /= it_scaled_hist->second;
            }
            hist.Ptr->Scale(scale);
            Scaled_Hists[hist.Ptr] = Global_MC_Scaling;
        }
    }
}

Long64_t hstack::Merge(TCollection* li)
{
    if(!li)
        return 0;
    TIter next(li);

    while(auto h = dynamic_cast<hstack*>(next())) {
        if(!IsCompatible(*h)) {
            LOG(ERROR) << "Skipping incompatible hstack:\n "
                       << *this << "\n"
                       << *h;
            continue;
        }
        // we add non-existing paths
        auto have_path = [] (const hists_t& hists, const string& path) {
            for(const auto& hist : hists)
                if(hist.Path == path)
                    return true;
            return false;
        };
        for(const auto& hist : h->hists) {
            if(!have_path(hists, hist.Path))
                hists.emplace_back(hist);
        }
    }

    return 0;
}

// cereal business

template<class Archive>
void save(Archive& archive,
          const TNamed& m)
{
    archive(string(m.GetName()), string(m.GetTitle()));
}

template<class Archive>
void load(Archive & archive,
          TNamed& m)
{
    string name;
    string title;
    archive(name, title);
    m.SetName(name.c_str());
    m.SetTitle(title.c_str());
}

namespace cereal
{
template <class Archive>
struct specialize<Archive, hstack, cereal::specialization::member_serialize> {};
}

void hstack::Streamer(TBuffer& R__b)
{
    stream_TBuffer buf(R__b);
    iostream inoutstream(addressof(buf));

    if (R__b.IsReading()) {
        cereal::BinaryInputArchive ar(inoutstream);
        ar(*this);
    }
    else {
        cereal::BinaryOutputArchive ar(inoutstream);
        ar(*this);
    }
}




