#include "plot/root_draw.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TObject.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "THStack.h"
#include "TVirtualPad.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>


using namespace ant;
using namespace std;

const padoption padoption::Legend = [] (TVirtualPad* p) {p->BuildLegend();};
const padoption padoption::LogX = [] (TVirtualPad* p) {p->SetLogx();};
const padoption padoption::LogY = [] (TVirtualPad* p) {p->SetLogy();};
const padoption padoption::LogZ = [] (TVirtualPad* p) {p->SetLogz();};

unsigned int canvas::num = 0;

const endcanvas ant::endc = endcanvas();
const endrow ant::endr = endrow();
const samepad_t ant::samepad = samepad_t();


canvas::canvas(const string& title) :
  name(), pads(), current_drawoption(), global_padoptions()
{
    CreateTCanvas(title);
}

canvas::~canvas()
{
}

TCanvas* canvas::CreateTCanvas(const string& title)
{
    name = std_ext::formatter() << "_canvas_" << setfill('0') << setw(3) << num++;
    return new TCanvas(name.c_str(), title.c_str());
}


TCanvas* canvas::FindTCanvas()
{
    TObject* o = gROOT->FindObjectAny(name.c_str());
    TCanvas* c = dynamic_cast<TCanvas*>(o);
    if(c)
        return c;
    else
        return CreateTCanvas();
}

void canvas::DrawObjs(TCanvas* c, unsigned cols, unsigned rows)
{
    c->Divide(cols,rows);
    int pad=1;
    unsigned ninrow =0;
    for(const pad_t& p : pads) {

        if(p.DrawableItems.empty()) {
            pad += (cols-ninrow);
            ninrow = 0;
        } else {

            TVirtualPad* vpad = c->cd(pad);
            // draw the objects
            for(const auto& item : p.DrawableItems) {
                item.Drawable->Draw(item.Option);
            }
            // set pad options
            for(const auto& option : p.PadOptions) {

                (*option)(vpad);
            }

            ++pad;
        }
        ninrow++;
    }
}

void canvas::cd()
{
    TCanvas* c = FindTCanvas();
    if(c) {
        c->cd();
    }
}


void canvas::AddDrawable(std::unique_ptr<root_drawable_traits> drawable)
{
    onetime_padoptions.insert(onetime_padoptions.end(),
                              global_padoptions.begin(),
                              global_padoptions.end());
    string drawoption = current_drawoption;
    if(!addobject || pads.empty()) {
        pads.emplace_back(onetime_padoptions);
    }
    else {
        drawoption += "same";
    }
    pads.back().DrawableItems.emplace_back(move(drawable), drawoption);

    onetime_padoptions.clear();
    addobject = false;
}

canvas& canvas::operator<<(root_drawable_traits& drawable)
{
    using container_t = drawable_container<root_drawable_traits*>;
    AddDrawable(std_ext::make_unique<container_t>(addressof(drawable)));
    return *this;
}

canvas& canvas::operator<<(TObject* hist)
{
    using container_t = drawable_container<TObject*>;
    AddDrawable(std_ext::make_unique<container_t>(hist));
    return *this;
}


canvas& canvas::operator<<(const endrow&)
{
    // just an empty pad indicates the end of the row
    pads.emplace_back();
    automode = false;
    return *this;
}

canvas& canvas::operator<<(const samepad_t&)
{
    addobject = true;
    return *this;
}

canvas& canvas::operator<<(const drawoption& c)
{
    current_drawoption = c.Option();
    return *this;
}

canvas&canvas::operator<<(const padoption& c)
{
    onetime_padoptions.emplace_back(std::addressof(c));
    return *this;
}

canvas& canvas::operator<<(const padoption::enable& c)
{

    global_padoptions.emplace_back(c.Modifier);
    return *this;
}

canvas& canvas::operator<<(const padoption::disable& c)
{
    global_padoptions.remove(c.Modifier);
    return *this;
}

canvas& canvas::operator<<(const endcanvas&)
{
    if(pads.empty()) {
        return *this;
    }

    TCanvas* c = FindTCanvas();

    if(c) {

        unsigned cols =0;
        unsigned rows =0;

        if(automode) {

            cols = ceil(sqrt(pads.size()));
            rows = ceil((double)pads.size()/(double)cols);

        } else {
            unsigned ccols=0;
            for(const auto& o : pads) {
                if(o.DrawableItems.empty()) {
                    cols = max(ccols,cols);
                    ccols=0;
                    rows++;
                } else {
                    ccols++;
                }
            }
        }
        DrawObjs(c,cols,rows);
    }

    return *this;
}

canvas& canvas::operator>>(const string& filename)
{
    TCanvas* c = FindTCanvas();
    if(c) {
        c->SaveAs(filename.c_str());
    }
    return *this;
}

//************ hstack* ***********


hstack::hstack(const string& name, const std::string& title,
               bool useIntelliLegend,
               bool ignoreEmptyHist):
    stack(new THStack(name.c_str(),title.c_str())),
    current_option(""),
    UseIntelliLegend(useIntelliLegend),
    IgnoreEmptyHist(ignoreEmptyHist)
{}

hstack::~hstack()
{}

hstack& hstack::operator<<(TH1* hist)
{
    hists.emplace_back(hist, current_option);
    xlabel = hist->GetXaxis()->GetTitle();
    ylabel = hist->GetYaxis()->GetTitle();
    return *this;
}

hstack& hstack::operator<<(const drawoption& c)
{
    current_option = c.Option();
    return *this;
}

void hstack::Draw(const string& option) const
{

    vector<string> orig_titles;
    if(UseIntelliLegend) {
        const auto& delim = ": ";
        vector<vector<string>> title_parts;
        map<string, unsigned>  token_counter;
        for(const auto& hist : hists) {
            const string title = hist.Ptr->GetTitle();
            orig_titles.emplace_back(title);
            title_parts.emplace_back(std_ext::tokenize_string(title, delim));
            for(const auto& token : title_parts.back())
                token_counter[token]++;
        }
        vector<vector<string>> unique_title_parts;
        for(size_t i=0; i< title_parts.size(); i++) {
            const auto& tokens = title_parts[i];
            unique_title_parts.emplace_back();
            auto& unique_tokens = unique_title_parts.back();
            for(const auto& token : tokens)
                if(token_counter[token] < hists.size())
                    unique_tokens.emplace_back(token);
            hists[i].Ptr->SetTitle(std_ext::concatenate_string(unique_tokens, delim).c_str());
        }
    }

    for(const auto& hist : hists) {
        if(IgnoreEmptyHist && hist.Ptr->GetEntries()==0)
            continue;
        stack->Add(hist.Ptr, hist.Option.c_str());
    }

    stack->Draw(option.c_str());

    auto xaxis = stack->GetXaxis();
    if(xaxis)
        xaxis->SetTitle(xlabel.c_str());

    auto yaxis = stack->GetYaxis();
    if(yaxis)
        yaxis->SetTitle(ylabel.c_str());

    if(UseIntelliLegend) {
        gPad->BuildLegend();
        for(size_t i=0;i<orig_titles.size();i++)
            hists[i].Ptr->SetTitle(orig_titles[i].c_str());
    }
}

const std::vector<Color_t> ColorPalette::Colors = {kRed, kGreen+1, kBlue, kYellow+1, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
