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
#include "TPaveText.h"

#include <list>

using namespace ant;
using namespace std;

hstack::options_t hstack::GlobalOptions;

double hstack::Global_MC_Scaling = std_ext::NaN;
map<TH1*, double> hstack::Scaled_Hists = {};

interval<double> hstack::GlobalYAxisRange = {std_ext::NaN, std_ext::NaN};
interval<interval<double>> hstack::GlobalLegendPosition = {
    {0.7,  0.63},
    {0.88, 0.88}
};


hstack::hstack(const string& name, const std::string& title) :
    TNamed(name.c_str(), title.c_str())
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
            ylabel == other.ylabel;
}

template<typename Att>
void AttCopy(const Att* src, Att* dest) {
    src->Copy(*dest);
}

bool isDataHist(const string& option)
{
    /// \todo find better way to detect which histograms
    /// should NOT be scaled or drawn with error bar in legend.
    /// For now, we assume that data is always drawn
    /// with error bars.
    return std_ext::contains(option, "E");
}



void hstack::buildIntelliTitle() const
{
    vector<string> title_parts = std_ext::tokenize_string(GetTitle(), ": ");
    if(title_parts.size()<3)
        return;

    const auto height = 0.03*(title_parts.size()-1);
    const auto x1 = GlobalLegendPosition.Start().Start();
    const auto y1 = GlobalLegendPosition.Start().Stop()-height;
    const auto x2 = GlobalLegendPosition.Stop().Start();
    const auto y2 = GlobalLegendPosition.Start().Stop();

    TPaveText *pt = new TPaveText(x1,y1,x2,y2, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetBorderSize(1);

    //pt->AddText(("Cuts: "+title_parts.back()).c_str());
    for(auto it = next(title_parts.begin()); it != prev(title_parts.end()); ++it)
        pt->AddText(it->c_str());
    pt->Draw();
}

hstack::ModOption_t unpackModOption_t(const std::string& hexstr) {
    stringstream buf_hex(hexstr);
    stringstream buf;
    unsigned int c; // for correct overload call of >>
    while(buf_hex >> hex >> c)
        buf << (unsigned char)c;
    cereal::BinaryInputArchive ar(buf);
    hstack::ModOption_t modoption;
    ar(modoption);
    return modoption;
}

std::string packModOption_t(const hstack::ModOption_t& modoption) {
    stringstream buf;
    cereal::BinaryOutputArchive ar(buf);
    ar(modoption);
    // convert to string hex to avoid null byte problems...
    stringstream buf_hex;
    // pay attention to correct conversion to hex with >>
    for(unsigned char c : buf.str()) {
        buf_hex << hex << (unsigned int)c << " ";
    }
    return buf_hex.str();
}

// this helper extracts the z number from the options,

struct tmphist_t {
    tmphist_t(TObject* obj, const char* packedoption) :
        Ptr(dynamic_cast<TH1*>(obj)),
        Option(unpackModOption_t(packedoption))
    {}
    TH1* Ptr;
    hstack::ModOption_t Option;
    double Entries;
    bool operator< (const tmphist_t& other) const {
        return Option.Z < other.Option.Z;
    }
};

void updateIntelliLegend(TLegend* legend, const list<tmphist_t>& hists)
{
    const auto& p = hstack::GlobalLegendPosition;
    legend->SetX1NDC(p.Start().Start());
    legend->SetY1NDC(p.Start().Stop());
    legend->SetX2NDC(p.Stop().Start());
    legend->SetY2NDC(p.Stop().Stop());

    legend->Clear();

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
    auto it_hist = hists.begin();
    for(size_t i=0; i< title_parts.size(); i++) {
        const tmphist_t& hist = *it_hist;

        const auto& tokens = title_parts[i];
        unique_title_parts.emplace_back();
        auto& unique_tokens = unique_title_parts.back();
        for(const auto& token : tokens)
            if(token_counter[token] < hists.size())
                unique_tokens.emplace_back(token);
        string unique_title = std_ext::concatenate_string(unique_tokens, delim);
        if(hstack::GlobalOptions.ShowEntriesInLegend)
            unique_title += std_ext::formatter() << " (" << hist.Entries << ")";

        auto entry = legend->AddEntry((TObject*)0, unique_title.c_str(), isDataHist(hist.Option.DrawOption) ? "lpfe" : "lpf");
        AttCopy<TAttLine>(hist.Ptr, entry);
        AttCopy<TAttFill>(hist.Ptr, entry);
        AttCopy<TAttMarker>(hist.Ptr, entry);
        if(hist.Option.BkgColor>=0) {
            // during legend building, modify fill color
            entry->SetFillColor(hist.Option.BkgColor);
            entry->SetFillStyle(1001);
            entry->SetLineWidth(2);
        }

        ++it_hist;
    }
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

        TLegend* legend = nullptr;

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

            // we need to iterate over links, since we need the option string
            // attached to each item
            if(fHists->IsEmpty())
                return;
            auto firstlink = fHists->FirstLink();
            list<tmphist_t> tmp_hists;
            do {
                tmp_hists.emplace_back(firstlink->GetObject(), firstlink->GetOption());
            }
            while((firstlink = firstlink->Next()));

            // save for later
            auto tmp_hists_backup = tmp_hists;

            // remove all empty hists if requested
            if(fHistogram) {
                const auto xfirst = fHistogram->GetXaxis()->GetFirst();
                const auto xlast = fHistogram->GetXaxis()->GetLast();
                auto it_hist = tmp_hists.begin();
                while(it_hist != tmp_hists.end()) {
                    TH1* h = it_hist->Ptr;
                    it_hist->Entries = h->GetDimension() > 1 ?
                                           h->GetEntries() : h->Integral(xfirst, xlast);
                    if(GlobalOptions.IgnoreEmptyHist && it_hist->Entries == 0)
                        it_hist = tmp_hists.erase(it_hist);
                    else
                        ++it_hist;
                }
            }

            // build the legend before re-ordering and after removal
            if(GlobalOptions.UseIntelliLegend) {
                if(!legend)
                    legend = new TLegend();
                updateIntelliLegend(legend, tmp_hists);
            }
            else {
                delete legend;
                legend = nullptr;
            }


            // reverse first
            std::reverse(tmp_hists.begin(), tmp_hists.end());

            // then sort according to Z (the sort is stable)
            tmp_hists.sort();


            fHists->Clear("nodelete");
            for(const auto& o : tmp_hists) {
                fHists->Add(o.Ptr, o.Option.DrawOption.c_str());
            }

            // let's hope that this does not throw an exception
            // before we unreverse the fHists
            THStack::Paint(chopt);

            if(legend)
                legend->Draw();

            // restore from backup
            fHists->Clear("nodelete");
            for(const auto& o : tmp_hists_backup) {
                fHists->Add(o.Ptr, packModOption_t(o.Option).c_str());
            }

        }

    };

    auto stack = new THStack_hack((string(GetName())+"_").c_str(),
                                  GlobalOptions.UseIntelliTitle ? "" : GetTitle());

    for(const auto& hist : hists)
        stack->Add(hist.Ptr, packModOption_t(hist.Option).c_str());

    UpdateMCScaling();

    // finally draw the stack
    string option_str(option);
    option_str += "nostack"; // always draw with nostack
    stack->Draw(option_str.c_str());

    // axis business
    auto xaxis = stack->GetXaxis();
    xaxis->SetTitle(xlabel.c_str());

    auto yaxis = stack->GetYaxis();
    yaxis->SetTitle(ylabel.c_str());

    if(GlobalOptions.UseIntelliTitle)
        buildIntelliTitle();
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
    gPad->Modified();
    gPad->Update();
}

void hstack::SetGlobalLegendPosition(double x1, double y1, double x2, double y2)
{
    GlobalLegendPosition = {{x1, y1}, {x2, y2}};
    gPad->Modified();
    gPad->Update();
}


void hstack::UseIntelliLegend(bool flag) { GlobalOptions.UseIntelliLegend = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetUseIntelliLegend() const { return GlobalOptions.UseIntelliLegend; }
void hstack::UseIntelliTitle(bool flag) { GlobalOptions.UseIntelliTitle = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetUseIntelliTitle() const { return GlobalOptions.UseIntelliTitle; }
void hstack::IgnoreEmptyHist(bool flag) { GlobalOptions.IgnoreEmptyHist = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetIgnoreEmptyHist() const { return GlobalOptions.IgnoreEmptyHist; }
void hstack::ShowEntriesInLegend(bool flag) { GlobalOptions.ShowEntriesInLegend = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetShowEntriesInLegend() const { return GlobalOptions.ShowEntriesInLegend; }

void hstack::UpdateMCScaling()
{
    if(isfinite(Global_MC_Scaling)) {
        for(hist_t& hist : hists) {
            if(isDataHist(hist.Option.DrawOption))
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




