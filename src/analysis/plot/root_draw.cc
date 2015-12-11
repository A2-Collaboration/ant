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

const padoption::map_options_t padoption::map_options =
{
    {padoption_t::Legend, [] (TVirtualPad* p) {p->BuildLegend();} },
    {padoption_t::LogX, [] (TVirtualPad* p) {p->SetLogx();} },
    {padoption_t::LogY, [] (TVirtualPad* p) {p->SetLogy();} },
    {padoption_t::LogZ, [] (TVirtualPad* p) {p->SetLogz();} },
};


unsigned int canvas::num = 0;

const endcanvas ant::endc;
const endrow ant::endr;
const samepad_t ant::samepad;


canvas::canvas(const string& title) :
  name(), pads(), current_drawoption(), current_padoptions()
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
                const auto& it = padoption::map_options.find(option);
                // silently ignore not implemented pad options
                if(it == padoption::map_options.end())
                    continue;
                it->second(vpad);
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
    if(!addobject || pads.empty()) {
        pads.emplace_back(current_padoptions);
    }
    pads.back().DrawableItems.emplace_back(move(drawable), current_drawoption);

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

canvas& canvas::operator<<(const padoption::set& c)
{
    const auto& o = c.Option();
    auto it = std::find(current_padoptions.begin(), current_padoptions.end(), o);
    if(it == current_padoptions.end()) {
        current_padoptions.emplace_back(o);
    }
    return *this;
}

canvas& canvas::operator<<(const padoption::unset& c)
{
    const auto& o = c.Option();
    current_padoptions.remove(o);
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


hstack::hstack(const string& name, const std::string& title):
    stack(new THStack(name.c_str(),title.c_str())),
    current_option("")
{}

hstack::~hstack()
{}

hstack& hstack::operator<<(TH1D* hist)
{
    stack->Add(hist, current_option.c_str());
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
    stack->Draw(option.c_str());

    auto xaxis = stack->GetXaxis();
    if(xaxis)
        xaxis->SetTitle(xlabel.c_str());

    auto yaxis = stack->GetYaxis();
    if(yaxis)
        yaxis->SetTitle(ylabel.c_str());
}

const std::vector<Color_t> ColorPalette::Colors = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
