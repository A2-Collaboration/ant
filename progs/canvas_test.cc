#include <iostream>

#include "TCanvas.h"
#include "TLine.h"
#include "TRint.h"
#include "TH1D.h"
#include "TLatex.h"

using namespace std;

class GUIIndicator {
public:
    virtual void SetPosition(double p) =0;
    virtual double GetPosition() const =0;
};

class IndicatorLine: public TLine, public GUIIndicator {
protected:
    class ILabel: public TLatex {
    public:
        using TLatex::TLatex;

        virtual void Delete(Option_t*) override {}
    };

    ILabel label;

    virtual void updateLabel()=0;

public:
    IndicatorLine(const std::string& caption=""):TLine(),
        label(fX1,fY1,caption.c_str())
    {
        SetBit(kCanDelete, false);
        SetVertical(kTRUE);
        SetLineColor(kBlue);
        SetLineWidth(3);
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
    virtual void SetX1(Double_t x1) override {fX1=x1;fX2=x1; updateLabel();}
    virtual void SetX2(Double_t x2) override {fX2=x2;fX1=x2; updateLabel();}

    // ignore all y movements
    virtual void SetY1(Double_t) override {updateLabel();}
    virtual void SetY2(Double_t) override {updateLabel();}

    // set y positions
    virtual void SetupY(Double_t y1, Double_t y2) {fY1 = y1; fY2=y2; updateLabel(); }

    virtual void SetPosition(double p) override { SetX1(p); }
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
    virtual void SetY1(Double_t y1) override {fY1=y1;fY2=y1; updateLabel();}
    virtual void SetY2(Double_t y2) override {fY2=y2;fY1=y2; updateLabel();}

    // ignore all x movements
    virtual void SetX1(Double_t) override {updateLabel();}
    virtual void SetX2(Double_t) override {updateLabel();}

    // set x positions
    virtual void SetupX(Double_t x1, Double_t x2) {fX1 = x1; fX2=x2; updateLabel(); }

    virtual void SetPosition(double p) override { SetY1(p); }
    virtual double GetPosition() const override { return GetY1(); }

};


class CalCanvas : public TCanvas {
public:
    CalCanvas(const std::string& name):
        TCanvas(name.c_str()) {
    }

    virtual void Draw(TH1* h) {
        this->cd();
        h->Draw();
        Update();
    }

    void SetMarker() {

        double x1,y1,x2,y2;
        GetRangeAxis(x1,y1,x2,y2);

        cout << x1 << " " << y1 << " / "  << x2 << " " << y2 << endl;


        VerticalIndicatorLine* tmp = new VerticalIndicatorLine("#pi^{0} pos");
        tmp->SetPosition((x1+x2)/2.0);
        tmp->SetupY(y1,y2);
        tmp->Draw();

        HorizontalIndicatorLine* tmph = new HorizontalIndicatorLine("#pi^{0} pos 2");
        tmph->SetPosition((y1+y2)/2.0);
        tmph->SetupX(x1,x2);
        tmph->Draw();

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

