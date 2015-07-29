#include <iostream>
#include <map>
#include <list>

#include "TCanvas.h"
#include "TLine.h"
#include "TRint.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TF1.h"

using namespace std;

class VirtualKnob {
public:
    std::string name = "";
    virtual double get() const =0;
    virtual void set(double v) =0;
};

class VirtualFunction {
public:
    virtual void Draw() =0;
    virtual std::list<VirtualKnob*> getKnobs() =0;
};

class FGaus: public VirtualFunction {
public:
    TF1* f = new TF1("","gaus",-10,10);

    class MyKnob: public VirtualKnob {
    public:
        TF1* f = nullptr;
        int p = 0;
        MyKnob(const std::string& n, TF1* func, int par):f(func),p(par) {
            name = n;
        }

        virtual double get() const override {
            return f->GetParameter(p);
        }

        virtual void set(double a) override {
            f->SetParameter(p,a);
        }

    };

    class MyWKnob: public VirtualKnob {
    public:
        TF1* f = nullptr;
        int p = 0;
        MyWKnob(const std::string& n, TF1* func):f(func) {
            name = n;
        }

        virtual double get() const override {
            return f->GetParameter(1) + f->GetParameter(2);
        }

        virtual void set(double a) override {
            f->SetParameter(2,a - f->GetParameter(1));
        }

    };

    MyKnob A = MyKnob("A",f,0);
    MyKnob x0 = MyKnob("x0",f,1);
    MyWKnob w = MyWKnob("w",f);
    std::list<VirtualKnob*> knobs = { &A, &x0, &w };

    FGaus() {
        f->SetParameter(0,10);
        f->SetParameter(1,1);
        f->SetParameter(2,3);
    }

    virtual void Draw() override {
        f->Draw("same");
    }
    virtual std::list<VirtualKnob*> getKnobs() {
        return knobs;
    }


};

class Fucntion: public VirtualFunction {
public:
    double values[2] = {};


    class MyKnobA: public VirtualKnob {
    public:
        double* v = nullptr;
        MyKnobA(const std::string& n, double* vp):v(vp) {
            name = n;
        }

        virtual double get() const override {
            return v[0];
        }

        virtual void set(double a) override {
            v[0] = a;
            v[1] = v[1]+1;
        }

    };

    class MyKnobB: public VirtualKnob {
    public:
        double* v = nullptr;
        MyKnobB(const std::string& n, double* vp):v(vp) {
            name = n;
        }

        virtual double get() const override {
            return v[1];
        }

        virtual void set(double a) override {
            v[1] = a;
        }

    };

    MyKnobA ka = MyKnobA("A", values);
    MyKnobB kb = MyKnobB("B", values);

    std::list<VirtualKnob*> knobs = {
        &ka, &kb
    };

    virtual std::list<VirtualKnob*> getKnobs() override { return knobs; }

    Fucntion() {}

    virtual void Draw() override {
        cout << "a: " << values[0] << " b: " << values[1] << endl;
    }

};

class update_notify_traits {
public:
    virtual void update_me() =0;
};

class GUIIndicator: public update_notify_traits {
public:
    virtual void SetPosition(double p) =0;
    virtual double GetPosition() const =0;
};



class IndicatorLine: public TLine, public GUIIndicator {
protected:
    VirtualKnob* knob = nullptr;

    class ILabel: public TLatex {
    public:
        using TLatex::TLatex;

        virtual void Delete(Option_t*) override {}
    };

    ILabel label;

    virtual void updateLabel()=0;

public:
    virtual void update_other() {
        if(update)
            update->update_me();
    }

    virtual void update_me() override {
        SetPosition(knob->get());
    }

    update_notify_traits* update = nullptr;
    IndicatorLine(VirtualKnob* k): TLine(),
        knob(k),
        label(fX1,fY1,knob->name.c_str())
    {
        SetBit(kCanDelete, false);
        SetVertical(kTRUE);
        SetLineColor(kBlue);
        SetLineWidth(3);
        fX1 = knob->get();
        fX2 = knob->get();
    }

    virtual void Delete(Option_t*) override {
    }

    virtual void Draw(Option_t *option="") override {
        updateLabel();
        TLine::Draw(option);
        label.Draw();
    }

};

class VerticalIndicatorLine: public IndicatorLine {
public:
protected:
    virtual void updateLabel() override {
            label.SetTextColor(GetLineColor());
            label.SetX(GetX1());
            label.SetY(GetY1());
            label.SetTextAngle(90);
    }

public:
    using IndicatorLine::IndicatorLine;
    virtual ~VerticalIndicatorLine() = default;

    // move both points the same way in x
    virtual void SetX1(Double_t x1) override {fX1=x1;fX2=x1; updateLabel(); knob->set(x1); update_other();}
    virtual void SetX2(Double_t x2) override {fX2=x2;fX1=x2; updateLabel(); knob->set(x2); update_other();}

    // ignore all y movements
    virtual void SetY1(Double_t) override {updateLabel();}
    virtual void SetY2(Double_t) override {updateLabel();}

    // set y positions
    virtual void SetupY(Double_t y1, Double_t y2) {fY1 = y1; fY2=y2; updateLabel(); }

    virtual void SetPosition(double p) override { fX1=p;fX2=p; }
    virtual double GetPosition() const override { return GetX1(); }

};

class HorizontalIndicatorLine: public IndicatorLine {
public:
protected:
    virtual void updateLabel() override {
            label.SetTextColor(GetLineColor());
            label.SetX(GetX1());
            label.SetY(GetY1());
            label.SetTextAngle(0);
    }

public:
    using IndicatorLine::IndicatorLine;

    virtual ~HorizontalIndicatorLine() = default;

    // move both points the same way in y
    virtual void SetY1(Double_t y1) override {fY1=y1;fY2=y1; updateLabel(); knob->set(y1);update_other();}
    virtual void SetY2(Double_t y2) override {fY2=y2;fY1=y2; updateLabel();knob->set(y2);update_other();}

    // ignore all x movements
    virtual void SetX1(Double_t) override {updateLabel();}
    virtual void SetX2(Double_t) override {updateLabel();}

    // set x positions
    virtual void SetupX(Double_t x1, Double_t x2) {fX1 = x1; fX2=x2; updateLabel(); }

    virtual void SetPosition(double p) override { fY1=p;fY2=p; }
    virtual double GetPosition() const override { return GetY1(); }

};


class CalCanvas : public TCanvas, public update_notify_traits {
public:

    FGaus f;
    std::list<GUIIndicator*> indicators;

    CalCanvas(const std::string& name):
        TCanvas(name.c_str()) {
    }

    virtual void Draw(TH1* h) {
        this->cd();
        h->Draw();
        Update();
    }

    virtual void update_me() override {
        f.Draw();
        for(auto& i : indicators) {
            i->update_me();
        }
        Update();
    }

    void SetMarker() {

        double x1,y1,x2,y2;
        GetRangeAxis(x1,y1,x2,y2);

        cout << x1 << " " << y1 << " / "  << x2 << " " << y2 << endl;

        for(auto& k : f.getKnobs()) {
            VerticalIndicatorLine* tmp = new VerticalIndicatorLine(k);
            tmp->update = this;
            tmp->SetupY(y1,y2);
            tmp->Draw();
            indicators.emplace_back(tmp);
        }
        f.Draw();
    }

    virtual void ShowGuidelines(TObject*, const Int_t, const char, const bool) override {}


};

int main(int argc, char** argv) {
    TRint app("omega",&argc,argv);

    CalCanvas c("test");
    TH1D* h = new TH1D("h","h",100,0,100);
    h->Fill(10,10);
    c.Draw(h);
    c.SetMarker();

    app.Run(kTRUE);

}

