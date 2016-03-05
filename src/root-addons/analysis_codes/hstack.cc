#include "hstack.h"



#include "base/std_ext/string.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
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
#include "TDirectory.h"
#include "TPaveText.h"

#include <list>
#include <iomanip>

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
    THStack(name.c_str(), title.c_str()),
    origtitle(title)
{
    gDirectory->Append(this);
}

hstack::hstack() : THStack() {}

hstack::~hstack() {}

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

bool hstack::wraphist_t::operator<(const hstack::wraphist_t& other) const {
    return Hist->Option.Z < other.Hist->Option.Z;
}

hstack& hstack::operator<<(TH1* hist)
{

    hists.emplace_back(hist, current_option);

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
            string(GetTitle()) == string(other.GetTitle());
}

template<typename Att>
void AttCopy(const Att* src, Att* dest) {
    src->Copy(*dest);
}

bool hstack::hist_t::isDataHist() const
{
    /// \todo find better way to detect which histograms
    /// should NOT be scaled or drawn with error bar in legend.
    /// For now, we assume that data is always drawn
    /// with error bars.
    return std_ext::contains(Option.DrawOption, "E");
}

void hstack::buildIntelliTitle()
{
    vector<string> title_parts = std_ext::tokenize_string(origtitle, ": ");
    if(title_parts.size()<3)
        return;

    const auto height = 0.03*(title_parts.size()-1);

    if(!intellititle) {
        const auto& p = GlobalLegendPosition;
        intellititle = std_ext::make_unique<TPaveText>(
                           p.Start().Start(),
                           p.Start().Stop()-height,
                           p.Stop().Start(),
                           p.Start().Stop(),
                           "NDC"
                           );
        intellititle->SetFillColor(kWhite);
        intellititle->SetBorderSize(1);
        //pt->AddText(("Cuts: "+title_parts.back()).c_str());
        for(auto it = next(title_parts.begin()); it != prev(title_parts.end()); ++it)
            intellititle->AddText(it->c_str());
    }

    if(intellilegend) {
        const auto x1 = intellilegend->GetX1NDC();
        const auto y1 = intellilegend->GetY1NDC()-height;
        const auto x2 = intellilegend->GetX2NDC();
        const auto y2 = intellilegend->GetY1NDC();
        intellititle->SetX1NDC(x1);
        intellititle->SetY1NDC(y1);
        intellititle->SetX2NDC(x2);
        intellititle->SetY2NDC(y2);
    }
}

void hstack::updateIntelliLegend(TLegend& legend, std::list<wraphist_t> wraphists)
{

    if(GlobalOptions.FixLegendPosition) {
        const auto& p = GlobalLegendPosition;
        legend.SetX1NDC(p.Start().Start());
        legend.SetY1NDC(p.Start().Stop());
        legend.SetX2NDC(p.Stop().Start());
        legend.SetY2NDC(p.Stop().Stop());
    }

    legend.Clear();

    const auto& delim = ": "; // this makes the tokens

    vector<vector<string>> title_parts;
    for(const auto& wraphist : wraphists) {
        const auto& hist = *wraphist.Hist;
        const auto& title = hist.Ptr->GetTitle();
        title_parts.emplace_back(std_ext::tokenize_string(title, delim));
    }

    auto it_hist = wraphists.begin();
    for(size_t i=0; i< title_parts.size(); i++) {
        const wraphist_t& wraphist = *it_hist;
        const hist_t& hist = *wraphist.Hist;

        vector<string> unique_tokens;
        const auto& tokens = title_parts[i];

        for(size_t j=0;j<tokens.size();j++) {
            const auto& token = tokens[j];
            bool equals_others = true;
            for(size_t k=0;k<title_parts.size();k++) {
                if(k==i)
                    continue;
                if(title_parts[k].size() != tokens.size())
                    continue;
                if(token != title_parts[k][j]) {
                    equals_others = false;
                    break;
                }
            }
            if(!equals_others)
                unique_tokens.push_back(token);
        }

        string unique_title = std_ext::concatenate_string(unique_tokens, delim);
        if(hstack::GlobalOptions.ShowEntriesInLegend)
            unique_title += std_ext::formatter() << " (" << setprecision(3) << wraphist.Entries << ")";

        auto entry = legend.AddEntry((TObject*)0, unique_title.c_str(), hist.isDataHist() ? "lpfe" : "lpf");
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

void hstack::Paint(const char* chopt)
{
    if(hists.empty())
        return;

    // ensure the histograms are there
    checkHists();

    UpdateMCScaling();

    if(GlobalOptions.UseIntelliTitle)
        SetTitle("");
    else
        SetTitle(origtitle.c_str());

    if(GlobalYAxisRange.IsSane()) {
        SetMinimum(GlobalYAxisRange.Start());
        SetMaximum(GlobalYAxisRange.Stop());
    }
    else {
        // resets it to default -1111, well, THStack knows what to do then
        SetMinimum();
        SetMaximum();
    }

    // make list of wrapped hists, in order to kick out elements and
    // count entries...
    list<wraphist_t> tmp_hists(hists.begin(), hists.end());

    // remove all empty hists if requested
    const auto basehist = fHistogram ? fHistogram : hists.front().Ptr;
    const auto xaxis = basehist->GetXaxis();
    auto it_hist = tmp_hists.begin();
    while(it_hist != tmp_hists.end()) {
        TH1* h = it_hist->Hist->Ptr;
        it_hist->Entries = h->GetDimension() > 1 ?
                               h->GetEntries() : h->Integral(xaxis->GetFirst(), xaxis->GetLast());
        if((GlobalOptions.IgnoreEmptyHist && it_hist->Entries < 1)
           || (std_ext::contains(h->GetTitle(), GlobalOptions.HistThresholdMatches)
               && it_hist->Entries < GlobalOptions.HistThreshold) )
            it_hist = tmp_hists.erase(it_hist);
        else
            ++it_hist;
    }

    // build the legend before re-ordering and after removal
    if(GlobalOptions.UseIntelliLegend) {
        if(!intellilegend) {
            const auto& p = GlobalLegendPosition;
            intellilegend = std_ext::make_unique<TLegend>(
                         p.Start().Start(),
                         p.Start().Stop(),
                         p.Stop().Start(),
                         p.Stop().Stop()
                         );
        }
        // sort by name if it's Bkg_
        // little hack to get legends stable...
        tmp_hists.sort([] (const wraphist_t& a, const wraphist_t& b) {
            if(a.Hist->Option.Z == b.Hist->Option.Z) {
                const string& title_a = a.Hist->Ptr->GetTitle();
                const string& title_b = b.Hist->Ptr->GetTitle();
                if(std_ext::contains(title_a,"Bkg_") && std_ext::contains(title_b,"Bkg_"))
                    return title_a < title_b;
                else
                    return false;
            }
            return false;
        });
        updateIntelliLegend(*intellilegend, tmp_hists);
    }
    else {
        intellilegend = nullptr;
    }

    if(GlobalOptions.UseIntelliTitle) {
        buildIntelliTitle();
    }
    else {
        intellititle = nullptr;
    }

    // reverse first
    std::reverse(tmp_hists.begin(), tmp_hists.end());

    // then sort according to Z
    tmp_hists.sort();


    if(fHists)
        fHists->Clear("nodelete");
    for(const auto& o : tmp_hists) {
        Add(o.Hist->Ptr, o.Hist->Option.DrawOption.c_str());
    }

    // let's hope that this does not throw an exception
    // before we unreverse the fHists
    string chopt_str(chopt);
    chopt_str += "nostack"; // always draw with nostack
    THStack::Paint(chopt_str.c_str());

    if(basehist != fHistogram) {
        const auto yaxis = basehist->GetYaxis();
        GetXaxis()->SetTitle(xaxis->GetTitle());
        GetYaxis()->SetTitle(yaxis->GetTitle());
        THStack::Paint(chopt_str.c_str());
    }

    if(intellilegend && !gPad->FindObject(intellilegend.get()))
       intellilegend->Draw();

    if(intellititle && !gPad->FindObject(intellititle.get()))
        intellititle->Draw();
}

void hstack::SetGlobalMCScaling(double scaling)
{
    if(isfinite(scaling) && scaling != 0) {
        Global_MC_Scaling = scaling;
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

void hstack::SetHistThreshold(double thresh, const char* matches)
{
    GlobalOptions.HistThreshold = thresh;
    GlobalOptions.HistThresholdMatches = matches;
    gPad->Modified();
    gPad->Update();
}

void hstack::FixLegendPosition(bool flag)
{
    TIter next(gPad->GetListOfPrimitives());
    while(TObject* obj = next()) {
        if(obj->IsA() == TLegend::Class()) {
            TLegend* legend = dynamic_cast<TLegend*>(obj);
            GlobalLegendPosition = {
                {legend->GetX1NDC(), legend->GetY1NDC()},
                {legend->GetX2NDC(), legend->GetY2NDC()}
            };
        }
    }
    GlobalOptions.FixLegendPosition = flag;
}
bool hstack::GetFixLegendPosition() const { return GlobalOptions.FixLegendPosition; }

void hstack::UseIntelliLegend(bool flag) { GlobalOptions.UseIntelliLegend = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetUseIntelliLegend() const { return GlobalOptions.UseIntelliLegend; }
void hstack::UseIntelliTitle(bool flag) { GlobalOptions.UseIntelliTitle = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetUseIntelliTitle() const { return GlobalOptions.UseIntelliTitle; }
void hstack::IgnoreEmptyHist(bool flag) { GlobalOptions.IgnoreEmptyHist = flag; gPad->Modified(); gPad->Update();}
bool hstack::GetIgnoreEmptyHist() const { return GlobalOptions.IgnoreEmptyHist; }
void hstack::ShowEntriesInLegend(bool flag) { GlobalOptions.ShowEntriesInLegend = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetShowEntriesInLegend() const { return GlobalOptions.ShowEntriesInLegend; }

void hstack::UpdateMCScaling()
{
    if(isfinite(Global_MC_Scaling)) {
        for(hist_t& hist : hists) {
            if(hist.isDataHist())
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

Long64_t hstack::Merge(TCollection* li, TFileMergeInfo*)
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






