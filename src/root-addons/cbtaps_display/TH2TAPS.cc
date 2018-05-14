#include "TH2TAPS.h"
#include "TMath.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include "TString.h"
#include "TH2DrawTool.h"
#include "TList.h"

using namespace std;
using namespace ant;

using vec = ant::matrixstack::Vector;


static TH2DrawTool::point_list MakeBaF2Shape();
static TH2DrawTool::point_list MakePbWO4Shape();

const Double_t b = 3.0;   // half distance of parallels of BaF2 in cm
const Double_t a = b*2/TMath::Sqrt(3); // edge length of BaF2 in cm


static const TH2DrawTool::point_list baf2_shape = MakeBaF2Shape();
static const TH2DrawTool::point_list pbwo4_shape = MakePbWO4Shape();

TH2TAPS::TH2TAPS(const string &name, const string &title): TH2Crystals(name,title)
{
    Build();
}

void TH2TAPS::FillElements(const TH2TAPS &h)
{
    TH2Crystals::FillElements(h);
}

void TH2TAPS::SetElements(const TH2TAPS &h)
{
    TH2Crystals::SetElements(h);
}

TH2DrawTool::point_list MakeBaF2Shape()
{
    matrixstack s;
    TH2DrawTool::point_list shape(7);

    const vec va(a,0);

    for(int i=0;i<6;++i) {
        shape.at(i) = s.Transform(va);
        s.Rotate(60.0*TMath::DegToRad());
    }
    shape.at(6) = shape.at(0);

    return shape;
}

TH2DrawTool::point_list MakePbWO4Shape()
{
    TH2DrawTool::point_list shape(5);

    shape.at(0) = vec(0,      0);
    shape.at(1) = vec(0,      b);
    shape.at(2) = vec(-a/2.0, b);
    shape.at(3) = vec(-a,     0);
    shape.at(4) = shape.at(0);

    return shape;
}

void TH2TAPS::DrawShape(TH2DrawTool& c, bool isBaF2) {
    if(isBaF2)
        c.Draw(baf2_shape);
    else
        c.Draw(pbwo4_shape);
    c.FinishShape();

    // FinishShape adds the bin, so we get it here again
    TH2PolyBin* bin_obj = (TH2PolyBin*) fBins->Last();
    TGraph* polygon = (TGraph*)bin_obj->GetPolygon();
    if(isBaF2)
        polygon->SetNameTitle("", Form("Element %d (BaF2)",fBins->GetEntries()-1));
    else
        polygon->SetNameTitle("", Form("Element %d (PbWO4)",fBins->GetEntries()-1));
}

void TH2TAPS::Build() {
    TH2DrawTool c(this);

    // invert the x-axis, to match the default coordinate system
    // this corresponds then to a "view to target" (not "from target")
    c.Scale(-1,1);

    // index s runs from sector 0..5
    // index r runs from ring 0..10 for each sector
    // index i runs from elements 0..r-1 for each ring

    const vec va(TMath::Sqrt(3)*b, b);

    for(size_t s=0;s<6;s++) {
        c.PushMatrix();
        //c.Translate(va);

        for(size_t r=0;r<11;r++) {
            c.Translate(va);
            c.PushMatrix();
            for(size_t i=0;i<r+1;i++) {

                if(r<2) {
                    // r<2 are PbWO4
                    c.PushMatrix();
                    c.Rotate(s*60*TMath::DegToRad()); // undo the sector rotation

                    // the four crystals can be obtained
                    // by inverting the x/y coordinates
                    DrawShape(c);
                    c.Scale(-1,1); // invert x => reflect at y-axis
                    DrawShape(c);
                    c.Scale(-1,-1); // add point reflection
                    DrawShape(c);
                    c.Scale(-1,1); // reflect at y-axis again
                    DrawShape(c);
                    c.PopMatrix();
                }
                else {
                    // r>=2 are BaF2
                    // but the most outer ring r==10
                    // has two crystal less than expected at sector dependent positions
                    const bool isCF = s==2 || s==5; // sector C and F are different
                    if(r<10 || ( !isCF && i>0 && i<r ) || (isCF && i>1) ) {
                        DrawShape(c, true);
                    }
                }

                c.Translate(vec(0,-2*b));
            }

            c.PopMatrix();
        }



        c.PopMatrix();
        c.Rotate(-60*TMath::DegToRad());
    }

    SetStats(kFALSE);
    SetXTitle("x [cm]");
    SetYTitle("y [cm]");

}
