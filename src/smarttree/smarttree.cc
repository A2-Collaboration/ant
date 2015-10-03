/**
 * @file smartttree.cc
 * @brief C++11 implementation ;-)
 *
 * @note When using TTree::Draw() histograms need to have a name when using bin specification:
 *       "expression>>name(bins,min,max)".
 *       if min == max, the range will be determined automatically by ROOT
 *        https://root.cern.ch/doc/master/classTTree.html#ac4016b174665a086fe16695aad3356e2
 */

#include "smarttree.h"

#include <map>
#include <list>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "base/interval.h"
#include "base/std_ext/string.h"

#include "TDirectory.h"
#include "TROOT.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TEventList.h"
#include "TRandom.h"


using namespace ant;
using namespace std;

class SmartTreeImpl: public SmartTree {
protected:
    mutable TTree* tree;

    map<std::string, interval<double>> range_cuts;
    static const interval<double> noRange;
    list<string> cuts;

    TCut cut;

    /**
     * @brief Test if a cut expression compliles
     * @param cut The cut expression
     * @return true=cut ok
     */
    bool TestCut(const string& cut);

    TCut buildCut() const;
    static void strAdd(ostream& stream, const string& branch, const interval<double>& i);

    list<string> canvases;

    bool autoUpdateEnabled = true;
    void AutoUpdate();

    static bool isSane(const string& expr) {
        if(std_ext::contains(expr,":")) {
            cout << "Error in expression " << expr << "\n";
            cout << "May not contain \":\". Use other Draw() command to specify more dimensions." << endl;
            return false;
        }
        return true;
    }

    void RemoveEventList() {
        auto list = tree->GetEventList();
        tree->SetEventList(nullptr);
        delete list;
    }


public:
    SmartTreeImpl(TTree *Tree, const string& name): SmartTree(name), tree(Tree) {
        gDirectory->Add(this);
    }

    virtual void Draw(const string &x, const int xbins=100) override;
    virtual void Draw(const string &x, const string& y, const int xbins=100, const int ybins=100) override;

    virtual void SetRange(const string& branch, double min, double max) override;
    virtual void SetRange(const string& expression, const interval<double> &range);
    virtual void RemoveRange(const string& branch) override;
    virtual void PrintCuts() const override;
    virtual void Update() override;
    virtual void RemoveAllCuts() override;

    virtual bool AddCut(const string &cut) override;
    virtual bool RemoveCut(const string &cut) override;

    virtual void SetAutoUpdate(bool update=true) override { autoUpdateEnabled = update; }
    virtual bool GetAutoUpdate() const override { return autoUpdateEnabled; }

    virtual void CloseAll() override;

    virtual const interval<double>& GetRange(const string& expression) const {
        const auto entry = range_cuts.find(expression);
        if(entry == range_cuts.end()) {
            return noRange;
        }
        return entry->second;
    }

    const TCut& GetCut() const { return cut; }

    TTree& GetTree() const { return *tree; }

    virtual ~SmartTreeImpl() {}
};

struct padstack {
    TVirtualPad* pad = nullptr;
    padstack(): pad(gPad) {}
    ~padstack() {
        if(pad)
            pad->cd();
    }
};


string getRandomString() {
    const unsigned r = floor(gRandom->Uniform(0,99999999));
    return "_rr" + to_string(r);
}

template <typename T>
T* GetObject(const string& name) {
    return dynamic_cast<T*>(gROOT->FindObject(name.c_str()));
}

class DrawCanvas: public TCanvas {
protected:
    string name;
    string smartree_name;
    static unsigned n;
    string drawoption;

    struct Viewport {
        interval<double> x = {0,0};
        interval<double> y = {0,0};
        Viewport() {}
    };

    Viewport getViewport() {
        Viewport p;
        GetRangeAxis(p.x.Start(),p.y.Start(),p.x.Stop(),p.y.Stop());
        return p;
    }

    virtual string buildDrawString(const SmartTreeImpl& smt) const =0;

public:
    DrawCanvas(const string& option): TCanvas(string("__c"+to_string(n)).c_str(),""),
        name("__h"+to_string(n++)), drawoption(option) {}
    virtual ~DrawCanvas() {}

    void SetSmartTreeName(const string& n) { smartree_name = n; }

    void Redraw(const SmartTreeImpl& smt) {

        padstack ps;
        this->cd();

        const auto ds = buildDrawString(smt);

        cout << ds << endl;

        smt.GetTree().Draw( ds.c_str(), smt.GetCut(), drawoption.c_str());

        Modified();
        Update();
    }

    void NotifyAxisChange(const string& expression, const interval<double>& range);

};

unsigned DrawCanvas::n = 0;

class DrawCanvas1D : public DrawCanvas {
private:
    string xExpr;
    int bins = 100;

    string buildDrawString(const SmartTreeImpl& smt) const override {
        const auto& xrange = smt.GetRange(xExpr);
        return std_ext::formatter() << xExpr << ">>" << name
                                    << "(" << bins  << "," << xrange.Start() << "," << xrange.Stop() << ")";
    }

    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;

public:
    DrawCanvas1D(const string& x, const int xbins=100): DrawCanvas(""), xExpr(x), bins(xbins) {}

    void SetBinsX(const int xbins) {
        bins = xbins;
    }

    virtual ~DrawCanvas1D() {}
};

class DrawCanvas2D: public DrawCanvas {
protected:
    string xExpr;
    string yExpr;
    int xbins;
    int ybins;

    string buildDrawString(const SmartTreeImpl& smt) const override {
        const auto xrange = smt.GetRange(xExpr);
        const auto yrange = smt.GetRange(yExpr);
        return std_ext::formatter() << yExpr << ":" << xExpr << ">>" << name
                                    << "(" << xbins  << "," << xrange.Start() << "," << xrange.Stop() << ","
                                           << ybins  << "," << yrange.Start() << "," << yrange.Stop() << ")";
    }

    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;

public:
    DrawCanvas2D(const string& x, const string& y, const int xBins=100, const int yBins=100):
        DrawCanvas("colz"), xExpr(x), yExpr(y), xbins(xBins), ybins(yBins) {}

    void SetBins(const int x, const int y) {
        xbins = x;
        ybins = y;
    }

    virtual ~DrawCanvas2D() {}
};


const interval<double> SmartTreeImpl::noRange = {-1,-1};

TCut SmartTreeImpl::buildCut() const
{
    stringstream s;
    bool first(true);
    for(const auto& cut : range_cuts) {

        if(!first) {
            s << "&& ";
        } else {
            first = false;
        }

        strAdd(s, cut.first, cut.second);

    }

    first = true;
    for(const auto& cut : cuts) {
        if(!first) {
            s << "&&";
        } else {
              first = false;
        }
        s << "(" << cut << ")";
    }

    return TCut(s.str().c_str());
}

void SmartTreeImpl::strAdd(ostream &stream, const string &branch, const interval<double> &i)
{
    stream <<  "(" << branch << ") >= " << i.Start() << " && "
           <<  "(" << branch << ") <= " << i.Stop()  << " ";
}

void SmartTreeImpl::AutoUpdate()
{
    if(autoUpdateEnabled)
        Update();
}

bool SmartTreeImpl::TestCut(const string &cut)
{
    const auto c = GetCut() + TCut(cut.c_str());
    cout << "Runngin gut " << c << endl;
    const auto evlist_name = getRandomString();
    const string drawstr = ">>"+evlist_name;

    const auto res = tree->Draw(drawstr.c_str(), c);

    if(res>=0) {
        auto list = GetObject<TEventList>(evlist_name);
        if(list) {
            auto old_list = tree->GetEventList();
            tree->SetEventList(list);
            delete old_list;
            const auto perc = list->GetN()*100.0 / tree->GetEntries();
            cout << "Event list set, " << list->GetN() << " entries (" << setprecision(4) << perc << "%)" << endl;
        }
    }
    return (res>=0);
}

void SmartTreeImpl::Draw(const string &x, const int xbins)
{
    if(!isSane(x)) {
        return;
    }
    auto canvas = new DrawCanvas1D(x,xbins);
    canvas->Redraw(*this);
    canvas->SetSmartTreeName(this->GetName());
    canvases.emplace_back(canvas->GetName());
}

void SmartTreeImpl::Draw(const string &x, const string& y, const int xbins, const int ybins)
{
    if(!isSane(x)) {
        return;
    }
    if(!isSane(y)) {
        return;
    }
    auto canvas = new DrawCanvas2D(x,y,xbins, ybins);
    canvas->Redraw(*this);
    canvas->SetSmartTreeName(this->GetName());
    canvases.emplace_back(canvas->GetName());
}

void SmartTreeImpl::SetRange(const string &branch, double min, double max)
{
    SetRange(branch, interval<double>(min,max));
}

void SmartTreeImpl::SetRange(const string &expression, const interval<double> &range)
{
    if((range.Start() == range.Stop()) || (range.Stop()<range.Start())) {
        RemoveRange(expression);
    } else {

        stringstream s;
        strAdd(s, expression, range);

        if(TestCut(s.str())) {
            range_cuts[expression] = range;
        }
    }

    AutoUpdate();
}

void SmartTreeImpl::RemoveRange(const string &branch)
{
    range_cuts.erase(branch);
    RemoveEventList();
    cut = buildCut();
    TestCut("");
    AutoUpdate();
}

void SmartTreeImpl::PrintCuts() const
{
    cout << "---- Active Cuts ----\n";
    for(const auto& cutentry : range_cuts) {
        cout << cutentry.first << ": " << cutentry.second << "\n";
    }
    for(const auto& cut : cuts) {
        cout << cut << "\n";
    }
    cout << "---------------------" << endl;
}

void SmartTreeImpl::Update()
{
    cut = buildCut();

    PrintCuts();

    auto i = canvases.cbegin();
    while(i!=canvases.cend()) {
        const string& canvas = *i;
        auto c = GetObject<DrawCanvas>(canvas);
        if(c) {
            c->Redraw(*this);
            ++i;
        } else {
            i = canvases.erase(i);
        }
    }
}

void SmartTreeImpl::RemoveAllCuts()
{
    range_cuts.clear();
    cuts.clear();
    RemoveEventList();
    AutoUpdate();
}

bool SmartTreeImpl::AddCut(const string &cut)
{
    if(TestCut(cut)) {
        cuts.emplace_back(cut);
        AutoUpdate();
        return true;
    }
    return false;
}

bool SmartTreeImpl::RemoveCut(const string &cut)
{
    const auto pos = find(cuts.begin(), cuts.end(), cut);

    if(pos != cuts.end()) {
        cuts.erase(pos);
        RemoveEventList();
        this->cut = buildCut();
        TestCut("");
        AutoUpdate();
        return true;
    }
    cout << "Cut " << cut << " not found" << endl;
    return false;
}

void SmartTreeImpl::CloseAll()
{
    for(const auto& canvas : canvases) {
        auto c = GetObject<DrawCanvas>(canvas);
        delete c;
    }

    canvases.clear();
}


SmartTree::~SmartTree()
{
}

SmartTree *SmartTree::Create(TTree *tree)
{
    return new SmartTreeImpl(tree, getRandomString());
}

void DrawCanvas::NotifyAxisChange(const string &expression, const interval<double> &range)
{
    auto st = GetObject<SmartTreeImpl>(smartree_name);
    if(!st) {
        cout << smartree_name << " not found" << endl;
        smartree_name = "";
        return;
    } else {
        st->SetRange(expression, range);
    }
}

void DrawCanvas1D::HandleInput(EEventType button, Int_t x, Int_t y)
{
    const Viewport old = getViewport();  //remove again and handle in draw

    TCanvas::HandleInput(button,x,y);

    if(button ==kButton1Up) {

        const Viewport newViewPort = getViewport();

        if(newViewPort.x != old.x) {
            cout << "Notify .." << endl;
            NotifyAxisChange(xExpr, newViewPort.x);
        }
    }
}

void DrawCanvas2D::HandleInput(EEventType button, Int_t x, Int_t y)
{
    const Viewport old = getViewport();  //remove again and handle in draw

    TCanvas::HandleInput(button,x,y);

    if(button ==kButton1Up) {

        const Viewport newViewPort = getViewport();

        if(newViewPort.x != old.x) {
            NotifyAxisChange(xExpr, newViewPort.x);
        }

        if(newViewPort.y != old.y) {
            NotifyAxisChange(yExpr, newViewPort.y);
        }
    }
}
