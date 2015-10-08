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



string getRandomString() {
    const unsigned r = floor(gRandom->Uniform(0,99999999));
    return "_rr" + to_string(r);
}

template <typename T>
T* GetObject(const string& name) {
    return dynamic_cast<T*>(gROOT->FindObject(name.c_str()));
}

struct padstack {
    TVirtualPad* pad = nullptr;

    padstack() noexcept: pad(gPad) {}

    ~padstack() {
        if(pad)
            pad->cd();
    }
};

struct TitleReset {
    string oldtitle;
    TCanvas* obj;
    TitleReset(TCanvas* n): oldtitle(n->GetTitle()), obj(n) {}
    ~TitleReset() { obj->SetTitle("SDSD"); }
};

struct CutStream : stringstream {

    using stringstream::stringstream;

    CutStream& add(const string& s) {
        if(tellp()!=0) {
            *this << " && ";
        }
        *this << s;
        return *this;
    }
};


class CutManager {
protected:
    list<string> cuts;
    string cut;

    string buildCut() const;
    void update();

    friend ostream& operator<<(ostream& stream, const CutManager& cutmgr);


public:
    CutManager() = default;
    virtual ~CutManager() = default;

    void Add(const string& c);
    void Remove(const string& c);
    void Pop();
    void RemoveAll();


    const string& Cut() const { return cut; }


};

ostream& operator<<(ostream& stream, const CutManager& cutmgr) {
    std::copy(cutmgr.cuts.cbegin(), cutmgr.cuts.cend(), std::ostream_iterator<string>(stream, "\n"));
    return stream;
}

string CutManager::buildCut() const
{
    CutStream newcut;
    for(const auto& c : cuts) {
        newcut.add(c);
    }

    return newcut.str();
}

void CutManager::update()
{
    cut = buildCut();
}

void CutManager::Add(const string &c)
{
    cuts.emplace_back(c);
    update();
}

void CutManager::Remove(const string &c)
{
    cuts.remove(c);
    update();
}

void CutManager::Pop()
{
    if(!cuts.empty()) {
        cuts.pop_back();
        update();
    }
}

void CutManager::RemoveAll()
{
    cuts.clear();
    update();
}

class RangeManager {
protected:
    map<string, interval<double>> ranges;

    static const interval<double> emptyRange;

public:
    RangeManager() = default;
    virtual ~RangeManager() = default;

    void Set(const string& expr, const interval<double> &range);
    void Remove(const string& expr);
    void RemoveAll();
    const interval<double>& Get(const string& expr) const;

};

const interval<double> RangeManager::emptyRange(-1,-1);

void RangeManager::Set(const string &expr, const interval<double> &range) {
    ranges[expr] = range;
}

void RangeManager::Remove(const string &expr) {
    ranges.erase(expr);
}

void RangeManager::RemoveAll() {
    ranges.clear();
}

const interval<double>& RangeManager::Get(const string &expr) const
{
    const auto entry = ranges.find(expr);
    if(entry == ranges.end()) {
        return emptyRange;
    }

    return entry->second;
}

class DrawCanvas;

class CanvasManager {
protected:
    list<string> canvas_names;

    DrawCanvas* GetCanvas(const string& name) {
        return GetObject<DrawCanvas>(name);
    }

public:
    CanvasManager() = default;
    virtual ~CanvasManager() = default;

    template <typename T, typename... Args_t>
    T* Add(Args_t&&... args) {
        auto canvas = new T(std::forward<Args_t>(args)...);
        auto name = getRandomString();
        canvas->SetName(name.c_str());
        canvas_names.emplace_back(move(name));
        return canvas;
    }

    template <typename Func>
    void Apply(Func f) {
        for(auto it=canvas_names.begin(); it!=canvas_names.end(); ) {
            auto canvas = GetCanvas(*it);
            if(!canvas) {
                it = canvas_names.erase(it);
            } else {
                f(canvas);
                ++it;
            }
        }
    }

    bool empty() const noexcept { return canvas_names.empty(); }

    void RemoveAll();
    void Remove(const string& name);

    DrawCanvas *AddByName(const string& name);
};

void CanvasManager::RemoveAll() {
    canvas_names.clear();
}

void CanvasManager::Remove(const string &name)
{
    canvas_names.remove(name);
}

DrawCanvas* CanvasManager::AddByName(const string &name)
{
    auto canvas = GetObject<DrawCanvas>(name);
    if(canvas)
        canvas_names.emplace_back(name);
    return canvas;
}


class EventListManager {
protected:
    struct state_t {
        TCut cut;
        Long64_t entries;
        state_t(const TCut& c="", const Long64_t e=-1) noexcept: cut(c), entries(e) {}
        ~state_t() = default;

        bool operator == (const state_t& other) const noexcept {
            return cut == other.cut && entries == other.entries;
        }

        bool operator != (const state_t& other) const noexcept {
            return !(*this == other);
        }
    };

    state_t last;

public:
    EventListManager() = default;
    ~EventListManager() = default;

    void Set(TTree& tree, const TCut &cut);
    void Reset(TTree& tree);
};

void EventListManager::Set(TTree &tree, const TCut &cut)
{
    state_t s(cut, tree.GetEntries());

    if(s != last) {
        cout << "Updating EventList..." << flush;
        auto old = tree.GetEventList();
        tree.SetEventList(nullptr);
        delete old;

        if(cut != "") {
            const auto name = getRandomString();
            const auto drawcmd = ">>"+name;
            tree.Draw(drawcmd.c_str(), cut);
            auto list = GetObject<TEventList>(name);
            if(list) {
                tree.SetEventList(list);
                last = s;
                cout << "Done" << endl;
            } else {
                cerr << "Error finding TEventList " << name << endl;
            }
        } else {
            cout <<" Done" << endl;
        }
    }
}

void EventListManager::Reset(TTree &tree) {
    tree.SetEventList(nullptr);
}


class SmartTreeImpl: public SmartTree {
protected:

    mutable TTree* tree;

    EventListManager evlmgr;

    CutManager cutmgr;
    RangeManager rangemgr;

    CanvasManager canvasmgr;


    bool TestCut(const string& cut);

        bool autoUpdateEnabled = true;
    void AutoUpdate();


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
    virtual void RemoveRange(const string& expr) override;
    virtual void PrintCuts() const override;
    virtual void Update() override;
    virtual void RemoveAllCuts() override;

    virtual bool Cut(const string &cut) override;
    virtual bool RemoveCut(const string &cut) override;
    virtual void UndoCut() override;
    virtual void Limit(Long64_t n_entries=-1) override;


    virtual TCut GetCut() const { return TCut(cutmgr.Cut().c_str()); }

    virtual void SetAutoUpdate(bool update=true) override { autoUpdateEnabled = update; }
    virtual bool GetAutoUpdate() const override { return autoUpdateEnabled; }

    virtual void CloseAll() override;

    virtual const interval<double>& GetRange(const string& expression) const {
        return rangemgr.Get(expression);
    }

    void UnlinkCanvas(const string& name);
    void RelinkCanvas(const string& name);

    TTree& GetTree() const { return *tree; }

    virtual ~SmartTreeImpl() {}
};





class DrawCanvas: public ant::SmartTreeCanvas {
protected:
    string histname;
    string smartree_name;
    static unsigned n;
    string drawoption;
    string localcut;

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
    virtual string buildTitle() const =0;
    virtual void postDraw() const =0;

    bool frozen = false;

public:
    DrawCanvas(const string& option): SmartTreeCanvas("__c"+to_string(n),""),
        histname("__h"+to_string(n++)), drawoption(option) {}
    virtual ~DrawCanvas() {}

    void SetSmartTreeName(const string& n) { smartree_name = n; SetTitle(buildTitle().c_str()); }

    void Redraw(const SmartTreeImpl& smt) {

        padstack ps;
        this->cd();

        const auto ds = buildDrawString(smt);

        cout << ds << endl;

        smt.GetTree().Draw( ds.c_str(), smt.GetCut() + TCut(localcut.c_str()), drawoption.c_str());

        postDraw();

        Modified();
        Update();
    }

    virtual void SetFrozen(bool f=true) {
        frozen = f;
        auto st = GetObject<SmartTreeImpl>(smartree_name);
        if(st) {
            if(frozen) {
                st->UnlinkCanvas(GetName());
            } else {
                st->RelinkCanvas(GetName());
            }
        }

    }

    virtual bool GetFrozen() const { return frozen; }

    void NotifyAxisChange(const string& expression, const interval<double>& range);

    virtual void SetCut(const char* cut) { localcut = cut; }
    virtual const char* GetCut() const { return localcut.c_str(); }

};

unsigned DrawCanvas::n = 0;

class DrawCanvas1D : public DrawCanvas {
private:
    string xExpr;
    int bins = 100;

    string buildDrawString(const SmartTreeImpl& smt) const override {
        const auto& xrange = smt.GetRange(xExpr);
        return std_ext::formatter() << xExpr << ">>" << histname
                                    << "(" << bins  << "," << xrange.Start() << "," << xrange.Stop() << ")";
    }

    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;

    virtual void postDraw() const override {
        TH1* h = GetObject<TH1>(histname);
        if(h) {
            h->SetXTitle(xExpr.c_str());
        }
    }

    virtual string buildTitle() const override {
        return xExpr + "  " + (GetFrozen() ? "[F]" : "   ");
    }

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
        return std_ext::formatter() << yExpr << ":" << xExpr << ">>" << histname
                                    << "(" << xbins  << "," << xrange.Start() << "," << xrange.Stop() << ","
                                           << ybins  << "," << yrange.Start() << "," << yrange.Stop() << ")";
    }

    virtual string buildTitle() const override {
        return xExpr + " vs " + yExpr + "  " + (GetFrozen() ? "[F]" : "   ");
    }

    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;

    virtual void postDraw() const override {
        TH2* h = GetObject<TH2>(histname);
        if(h) {
            h->SetXTitle(xExpr.c_str());
            h->SetYTitle(yExpr.c_str());
        }
    }

public:
    DrawCanvas2D(const string& x, const string& y, const int xBins=100, const int yBins=100):
        DrawCanvas("colz"), xExpr(x), yExpr(y), xbins(xBins), ybins(yBins) {}

    void SetBins(const int x, const int y) {
        xbins = x;
        ybins = y;
    }

    virtual ~DrawCanvas2D() {}
};

void SmartTreeImpl::AutoUpdate()
{
    if(autoUpdateEnabled) {

        Update();
    }
}

bool SmartTreeImpl::TestCut(const string &cut)
{
    const auto c = GetCut() + TCut(cut.c_str());

    const auto res = tree->Draw("",c,"",0);

    return (res>=0);
}

void SmartTreeImpl::Draw(const string &x, const int xbins)
{
    auto canvas = canvasmgr.Add<DrawCanvas1D>(x, xbins);
    canvas->SetSmartTreeName(GetName());
    if(autoUpdateEnabled) {
        evlmgr.Set(*tree, GetCut());
        canvas->Redraw(*this);
    }
}

void SmartTreeImpl::Draw(const string &x, const string& y, const int xbins, const int ybins)
{
    auto canvas = canvasmgr.Add<DrawCanvas2D>(x, y, xbins, ybins);
    canvas->SetSmartTreeName(GetName());
    if(autoUpdateEnabled) {
        evlmgr.Set(*tree,GetCut());
        canvas->Redraw(*this);
    }
}

void SmartTreeImpl::SetRange(const string &branch, double min, double max)
{
    SetRange(branch, interval<double>(min,max));
}

void SmartTreeImpl::SetRange(const string &expression, const interval<double> &range)
{
    if((range.Start() == range.Stop()) || (range.Stop()<range.Start())) {
        rangemgr.Remove(expression);
    } else {
        //TODO: TestCut
        rangemgr.Set(expression,range);
    }

    AutoUpdate();
}

void SmartTreeImpl::RemoveRange(const string &expr)
{
    rangemgr.Remove(expr);
    AutoUpdate();
}

void SmartTreeImpl::PrintCuts() const
{
    cout << cutmgr << endl;
}

void SmartTreeImpl::Update()
{
    PrintCuts();

    if(!canvasmgr.empty()) {
        evlmgr.Set(*tree, GetCut());



        canvasmgr.Apply([this] (DrawCanvas* c) {
            c->Redraw(*this);
        });
    }

}

void SmartTreeImpl::RemoveAllCuts()
{
    cutmgr.RemoveAll();
    AutoUpdate();
}

bool SmartTreeImpl::Cut(const string &cut)
{
    if(TestCut(cut)) {
        cutmgr.Add(cut);
        AutoUpdate();
        return true;
    }
    return false;
}

bool SmartTreeImpl::RemoveCut(const string &cut)
{
    //TODO: return something useful
    cutmgr.Remove(cut);
    return true;
}

void SmartTreeImpl::UndoCut()
{
    cutmgr.Pop();
    AutoUpdate();
}

void SmartTreeImpl::Limit(Long64_t n_entries)
{
    if(n_entries > tree->GetEntries()) {
        n_entries = tree->GetEntries();
    }

    if(n_entries < -1) {
        n_entries = -1;
    }

    tree->SetEntries(n_entries);
}

void SmartTreeImpl::CloseAll()
{
    canvasmgr.Apply([] (DrawCanvas* c) { delete c;});
    canvasmgr.RemoveAll();
}

void SmartTreeImpl::UnlinkCanvas(const string &name)
{
    canvasmgr.Remove(name);
}

void SmartTreeImpl::RelinkCanvas(const string &name)
{
    auto canvas = canvasmgr.AddByName(name);
    if(autoUpdateEnabled && canvas) {
        canvas->Redraw(*this);
    }
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
    if(smartree_name.empty())
        return;

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

SmartTreeCanvas::SmartTreeCanvas(const string &name, const string &title):
    TCanvas(name.c_str(), title.c_str())
{}

SmartTreeCanvas::~SmartTreeCanvas()
{

}
