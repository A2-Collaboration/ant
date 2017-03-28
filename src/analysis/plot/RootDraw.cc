#include "plot/RootDraw.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"

#include "TH1.h"
#include "TObject.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "THStack.h"
#include "TVirtualPad.h"
#include "TTree.h"

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
const padoption padoption::MakeSquare = [] (TVirtualPad* p) {
    auto min_size = std::min(p->GetWh(), p->GetWw());
    p->SetCanvasSize(min_size, min_size);
};

padoption::SetFillColor::SetFillColor(Color_t col) :
    padmodifier_t([col] (TVirtualPad* p) {p->SetFillColor(col);})
{}


unsigned int canvas::num = 0;

const endcanvas ant::endc = endcanvas();
const endrow ant::endr = endrow();
const samepad_t ant::samepad = samepad_t();


canvas::canvas(const string& title_) :
  name(std_ext::formatter() << "_canvas_" << setfill('0') << setw(3) << num++),
  title(title_),
  pads(), current_drawoption(), global_padoptions()
{
}

canvas::~canvas()
{
}

TCanvas* canvas::FindTCanvas() const
{
    TObject* o = gROOT->FindObjectAny(name.c_str());
    TCanvas* c = dynamic_cast<TCanvas*>(o);
    if(c)
        return c;
    else
        return new TCanvas(name.c_str(), title.c_str());
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
                option(vpad);
            }

            ++pad;
            ++ninrow;
        }
    }
}

void canvas::cd()
{
    TCanvas* c = FindTCanvas();
    if(c) {
        c->cd();
    }
}


void canvas::AddDrawable(std::shared_ptr<root_drawable_traits> drawable)
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

canvas& canvas::operator<<(TObject* hist)
{
    struct TObject_wrapper : root_drawable_traits {
        TObject* h;
        TObject_wrapper(TObject* h_) : h(h_) {}
        virtual ~TObject_wrapper() = default;
        virtual void Draw(const std::string& option) const {
            h->Draw(option.c_str());
        }
    };
    AddDrawable(std::make_shared<TObject_wrapper>(hist));
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

canvas&canvas::operator<<(const padmodifier_t& c)
{
    onetime_padoptions.emplace_back(c);
    return *this;
}

canvas& canvas::operator<<(const padoption::enable& c)
{
    global_padoptions.emplace_back(c.Modifier);
    return *this;
}

canvas& canvas::operator<<(const padoption::disable& c)
{
    // checking std::function for equality is a bit difficult
    global_padoptions.remove_if([c] (const padmodifier_t& mod) {
        typedef void (fnType)(TVirtualPad*);
        return mod.template target<fnType>() == c.Modifier.template target<fnType>();
    });
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
                // endr statements have no drawable items
                if(o.DrawableItems.empty()) {
                    cols = max(ccols,cols);
                    ccols=0;
                    rows++;
                } else {
                    ccols++;
                }
            }

            // if last statement wasn't endr, then add another row
            if(!pads.back().DrawableItems.empty())
                rows++;
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

const std::vector<Color_t> ColorPalette::Colors = {kRed, kGreen+1, kBlue, kYellow+1, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};

unsigned TTree_drawable::nInstances = 0;

TTree_drawable::TTree_drawable(TTree* tree, const string& formula, const string& cut) :
    Tree(tree), Cut(cut), Xlabel(""), Ylabel(""), autolabels(true)
{
    if(formula.empty()) {
        LOG(WARNING) << "Provided empty formula to TTree_drawable, that's strange...";
        Name = "";
        Title = "";
        Formula = "";
        nInstances++;
        return;
    }

    // Sanitize the formula provided
    const auto pos_shift_op = formula.find(">>");
    // search after >> operator for opening parenthesis, which indicates binning statement
    const auto pos_opening_parenthesis = formula.find("(", pos_shift_op);
    // between the >> operator and the binning, there should be the histname
    // use string_sanitize to get rid of possible leading/trailing whitespace
    const auto pos_behind_shift_op = pos_shift_op+2;
    const auto histname = std_ext::string_sanitize(formula.substr(pos_behind_shift_op, pos_opening_parenthesis-pos_behind_shift_op));
    // do not touch formulas with existing histnames
    if(histname != "") {
        Title = histname;
        Name = histname;
        Formula = formula;
        nInstances++;
        return;
    }

    // Last case, hist title is not specified, create generic
    Name = std_ext::formatter() << "htemp" << nInstances;
    const auto formula_clean = formula.substr(0, pos_shift_op);
    const auto binning_statement = pos_opening_parenthesis == string::npos ? "" : formula.substr(pos_opening_parenthesis);
    Formula = std_ext::formatter() << formula_clean << ">>" << Name << binning_statement;

    nInstances++;
    return;

}

TTree_drawable::TTree_drawable(TTree* tree, const string& varx, const string& cut, const string &title, const string &xlabel, const string &ylabel, const BinSettings &binsx, const std::string &name) :
    Tree(tree), Cut(cut), Title(title), Name(name), Xlabel(xlabel), Ylabel(ylabel), autolabels(false)
{
    if(name == "")
    {
        LOG(WARNING) << "Provided empty title to TTree_drawable, that's strange... overwriting to htemp" << nInstances;
        Title = std_ext::formatter() << "htemp" << nInstances;
    }

    if (varx.find(">>") != std::string::npos) LOG(ERROR) << "Pipe (>>) included in TTree_drawable : " << varx;

    Formula = std_ext::formatter() << varx
                                   << " >> "
                                   << Name
                                   << "(" << binsx.Bins() << ","  << binsx.Start() << "," << binsx.Stop() << ")";
    nInstances++;
    return;

}

TTree_drawable::TTree_drawable(TTree* tree, const string& varx, const string& vary, const string& cut, const string &title, const string &xlabel, const string &ylabel, const BinSettings &binsx, const BinSettings &binsy, const std::string &name) :
    Tree(tree), Cut(cut), Title(title), Name(name), Xlabel(xlabel), Ylabel(ylabel), autolabels(false)
{
    if(name == "")
    {
        LOG(WARNING) << "Provided empty title to TTree_drawable, that's strange... overwriting to htemp" << nInstances;
        Title = std_ext::formatter() << "htemp" << nInstances;
    }

    if (varx.find(">>") != std::string::npos) LOG(ERROR) << "Pipe (>>) included in TTree_drawable : " << varx;
    if (vary.find(">>") != std::string::npos) LOG(ERROR) << "Pipe (>>) included in TTree_drawable : " << vary;


    Formula = std_ext::formatter() << vary << ":" << varx
                                   << " >> "
                                   << Name
                                   << "(" << binsx.Bins() << ","  << binsx.Start() << "," << binsx.Stop() << ","
                                          << binsy.Bins() << ","  << binsy.Start() << "," << binsy.Stop() << ")";
    nInstances++;
    return;

}

void TTree_drawable::Draw(const string& option) const
{
    Tree->Draw(Formula.c_str(), Cut.c_str(), option.c_str());
    auto h = (TH1*)gDirectory->Get(Name.c_str());

    if((h) && !(autolabels))
    {
        h->SetTitle(Name.c_str());
        h->GetXaxis()->SetTitle(Xlabel.c_str());
        h->GetYaxis()->SetTitle(Ylabel.c_str());
    }

}
