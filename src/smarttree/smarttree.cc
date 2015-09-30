/**
 * @file smartttree.cc
 * @brief C++11 implementation ;-)
 */

#include "smarttree.h"

#include <map>
#include <list>
#include <iostream>
#include <sstream>

#include "base/interval.h"

#include "TDirectory.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TAxis.h"
#include "TCanvas.h"

using namespace ant;
using namespace std;

struct padstack {
    TVirtualPad* pad = nullptr;
    padstack(): pad(gPad) {}
    ~padstack() {
        if(pad)
            pad->cd();
    }
};

class DrawCanvas: public TCanvas {
protected:
    string drawstring;
    string drawoption;

public:
    DrawCanvas(const string& dstring, const string& option): drawstring(dstring), drawoption(option) {}
    virtual ~DrawCanvas() {}

    void Redraw(TTree& tree, const TCut& cut) {
        padstack ps;
        this->cd();
        tree.Draw(drawstring.c_str(), cut, drawoption.c_str());
    }
};

class DrawCanvas1D : public DrawCanvas {
public:
    DrawCanvas1D(const string& x): DrawCanvas(x,"") {}
    virtual ~DrawCanvas1D() {}
};

class DrawCanvas2D: public DrawCanvas {
public:
    DrawCanvas2D(const string& x, const string& y):
        DrawCanvas(y+":"+x,"colz") {}
    virtual ~DrawCanvas2D() {}
};

class SmartTreeImpl: public SmartTree {
protected:
    TTree* tree;

    map<std::string, interval<double>> range_cuts;
    list<string> cuts;

    /**
     * @brief Test if a cut expression compliles
     * @param cut The cut expression
     * @return true=cut ok
     */
    bool TestCut(const string& cut);

    TCut buildCut() const;
    static void strAdd(ostream& stream, const string& branch, const interval<double>& i);

    list<DrawCanvas*> canvases;

    bool autoUpdateEnabled = true;
    void AutoUpdate();



public:
    SmartTreeImpl(TTree *Tree): tree(Tree) {}

    virtual void Draw(const string &x) override;
    virtual void Draw(const string &x, const string& y) override;

    virtual void SetRange(const string& branch, double min, double max) override;
    virtual void RemoveRange(const string& branch) override;
    virtual void PrintCuts() const override;
    virtual void Update() override;
    virtual void RemoveAllCuts() override;

    virtual bool AddCut(const string &cut) override;
    virtual bool RemoveCut(const string &cut) override;

    virtual void SetAutoUpdate(bool update=true) override { autoUpdateEnabled = update; }
    virtual bool GetAutoUpdate() const override { return autoUpdateEnabled; }

    virtual ~SmartTreeImpl() {}
};


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

    for(const auto& cut : cuts) {
        s << "(" << cut << ") ";
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
    const auto c = buildCut() + TCut(cut.c_str());
    const auto res = tree->Draw("", c, "", 0);
    return (res>=0);
}

void SmartTreeImpl::Draw(const string &x)
{
    auto canvas = new DrawCanvas1D(x);
    canvas->Redraw(*tree, buildCut());
    canvases.emplace_back(canvas);
}

void SmartTreeImpl::Draw(const string &x, const string& y)
{
    auto canvas = new DrawCanvas2D(x,y);
    canvas->Redraw(*tree, buildCut());
    canvases.emplace_back(canvas);
}

void SmartTreeImpl::SetRange(const string &branch, double min, double max)
{

    if((min==max) || (min>max)) {
        RemoveRange(branch);
    } else {

        const interval<double> range(min,max);

        stringstream s;
        strAdd(s, branch, range);

        if(TestCut(s.str())) {
            range_cuts[branch] = range;
        }
    }

    AutoUpdate();
}

void SmartTreeImpl::RemoveRange(const string &branch)
{
    range_cuts.erase(branch);
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
    const auto cut = buildCut();

    PrintCuts();

    for(auto canvas : canvases ) {
        canvas->Redraw(*tree, cut);
    }
}

void SmartTreeImpl::RemoveAllCuts()
{
    range_cuts.clear();
    cuts.clear();
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
        AutoUpdate();
        return true;
    }
    cout << "Cut " << cut << " not found" << endl;
    return false;
}


SmartTree::~SmartTree()
{
}

SmartTree *SmartTree::Create(TTree *tree)
{
    return new SmartTreeImpl(tree);
}
