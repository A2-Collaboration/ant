#include "TH2CB.h"
#include "TText.h"
#include "Rtypes.h"
#include "TGraph.h"
#include "TH2DrawTool.h"
#include "TList.h"
#include "base/matrixstack.h"

#include "cb_numbering.h"

#include <iostream>

using namespace std;
using namespace ant;

using vec = ant::matrixstack::Vector;

const TH2DrawTool::point_list MakeTriangle();
const std::set<Int_t> MakeListOfBinsInholes();
const std::vector<Int_t> Make_binmap();
//void MakeLevel(TH2DrawTool &c, const UInt_t n, set<Int_t>::const_iterator &nexthole, Int_t& vbins);


// to be called in this order. List of holes has to be build before the bin map can be created.
static const TH2DrawTool::point_list shape = MakeTriangle();    // shape of a crystal
const std::set<Int_t> bins_in_holes = MakeListOfBinsInholes();
const std::vector<Int_t> binmap = Make_binmap();




void TH2CB::Build()
{
    TH2DrawTool tool(this);

    const vec& a(shape.at(1));
    const vec& b(shape.at(2));



    std::set<Int_t>::const_iterator nexthole = bins_in_holes.begin();
    Int_t vbins=0;
    tool.PushMatrix();
    tool.PushMatrix();
    for(int i=0;i<4;++i) {

        for( int j=0; j<2; ++j) {
            MakeLevel(tool,1,nexthole,vbins);
            tool.Translate(a);
            tool.Scale(-1,-1);
            MakeLevel(tool,1,nexthole,vbins);
            tool.Translate(b);
            tool.Scale(-1,-1);
        }
        tool.Translate(2.0*b-a);
    }
    tool.PopMatrix();

    tool.Translate(-1.0*a);

    for( int j=0; j<2; ++j) {
        MakeLevel(tool,1,nexthole,vbins);
        tool.Translate(a);
        tool.Scale(-1,-1);
        MakeLevel(tool,1,nexthole,vbins);
        tool.Translate(b);
        tool.Scale(-1,-1);
    }



    SetStats(kFALSE);
    GetXaxis()->SetTickLength(0);
    GetXaxis()->SetLabelSize(0);
    GetYaxis()->SetTickLength(0);
    GetYaxis()->SetLabelSize(0);

    TText* inlabel = new TText(3.24,-0.43,"Beam In");
    inlabel->SetTextSize(0.03);
    inlabel->SetTextAlign(22);  // middle, middle
    GetListOfFunctions()->Add(inlabel);

    TText* outlabel = new TText(0.73,-0.43,"Beam out");
    outlabel->SetTextSize(0.03);
    outlabel->SetTextAlign(22);  // middle, middle
    GetListOfFunctions()->Add(outlabel);

    for( int major=1;major<=20;++major){
        for( int minor=1;minor<=4;++minor) {
            for( int crystal=1;crystal<=9;++crystal) {

                const Int_t bin = GetBinOfMMC(major, minor, crystal);

                if( bin >0) {

                    TH2PolyBin* bin_obj = ((TH2PolyBin*) fBins->At(bin-1));
                    TGraph* graph = (TGraph*) bin_obj->GetPolygon();

                    const Int_t vbin = GetVBinOfMMC(major,minor,crystal);
                    if( vbin >=0 ) {
                        Int_t element = GetElementOfCrystal(vbin-1);

                        graph->SetNameTitle("", Form("Element %d (%d/%d/%d)",element, major, minor, crystal));
                    }
                }
            }

        }
    }

    tool.PopMatrix();

    // hacked in Glue pads for Xmas ball
    if(draw_glue_pads) {
        const auto ma = (TH2DrawTool::Vector)a*(1/6.0);
        const auto mb = (TH2DrawTool::Vector)b*(1/6.0);
        const auto mc = ma-mb;

        TH2DrawTool::point_list glue(4);
        glue.at(0) = TH2DrawTool::Vector(0.0, 0.0);
        glue.at(1) = ma;
        glue.at(2) = glue.at(1)+(TH2DrawTool::Vector)mc*4.0;
        glue.at(3) = glue.at(2)+(TH2DrawTool::Vector)(mc-mb);

        tool.PushMatrix();
            tool.Translate(b);
            tool.Translate((TH2DrawTool::Vector)a*-1);
            tool.PushMatrix();
            for(int i=0;i<5;++i) {
                tool.Draw(glue);
                tool.FinishShape();
                if(i!=4)tool.Translate(a);
            }
            tool.Translate((TH2DrawTool::Vector)mc*6.0);
            tool.Draw(glue);
            tool.FinishShape();
            tool.Translate(a);
            tool.PopMatrix();

        tool.PopMatrix();

        tool.PushMatrix();
            tool.Translate((TH2DrawTool::Vector)mc*6.0);
            tool.Translate((TH2DrawTool::Vector)b*-1);
            tool.Scale(1,-1);
            tool.PushMatrix();
            for(int i=0;i<5;++i) {
                tool.Draw(glue);
                tool.FinishShape();
                tool.Translate(a);
            }
            tool.PopMatrix();

        tool.PopMatrix();
    }

}

void TH2CB::MakeLevel(TH2DrawTool& c, const UInt_t n, set<Int_t>::const_iterator &nexthole, Int_t& vbins)
{

    const vec& a(shape.at(1));
    const vec& b(shape.at(2));

    if(n>=4) {

        ++vbins;

        if( (nexthole == bins_in_holes.end()) || (*nexthole!=vbins)) {
            c.Draw(shape);
            c.FinishShape();
        } else
            if ( nexthole != bins_in_holes.end() ) {
                   ++nexthole;
        }

    } else {
        c.PushMatrix();
        c.Scale(1.0/n,1.0/n);
       for( UInt_t row=0; row<n; ++row){
            const UInt_t triganles = 2*row +1;
            c.PushMatrix();
            c.Translate( (n-row-1.0)* b );
            for( UInt_t t=0; t<triganles; ++t) {
                if(t%2==0) {
                    MakeLevel(c,n+1,nexthole,vbins);
                    c.Translate(a);}
                else {
                    c.Translate(b);
                    c.Scale(-1,-1);
                    MakeLevel(c,n+1,nexthole,vbins);
                    c.Translate(b);
                    c.Scale(-1,-1);
                }
            }
            c.PopMatrix();
        }
        c.PopMatrix();
    }

}



TH2CB::TH2CB(const string &name, const string &title, bool glue_pads): TH2Crystals(name,title),
  draw_glue_pads(glue_pads)
{
    Build();
}

Int_t TH2CB::GetBinOfMMC(const UChar_t major, const UChar_t minor, const UChar_t crystal)
{
    const Int_t vbin = GetVBinOfMMC(major,minor,crystal);
    return binmap.at(vbin);
}

Int_t TH2CB::GetVBinOfMMC(const UChar_t major, const UChar_t minor, const UChar_t crystal)
{
    if (    ( major<1 || major> 20)
         || ( minor<1 || minor>4 )
         || ( crystal<1 || crystal>9 )  ) {
        return -5;  // invalid element specifications mapped to the "sea"
    }

    return (major-1)*4*9 + (minor-1)*9 + crystal;
}

Int_t TH2CB::GetBinOfVBin(const Int_t vbin)
{
    return binmap.at(vbin);
}

void TH2CB::FillCrystalNumbers()
{
    for(Int_t i=1; i<=720; ++i ) {
        const Int_t bin = GetBinOfVBin(i);
        if(bin>0)
            SetBinContent(bin,i-1);
    }
}

void TH2CB::FillMMCNumbers()
{
    for( int major=1;major<=20;++major){
        for( int minor=1;minor<=4;++minor) {
            for( int crystal=1;crystal<=9;++crystal) {
                const Int_t number = major*100 + minor*10+crystal;
                const Int_t bin = GetBinOfMMC(major,minor,crystal);
                if(bin>0)
                    SetBinContent(bin,number);
            }
        }
    }
}

void TH2CB::FillElementNumbers()
{
    for(Int_t i=0; i<720; ++i ) {
        SetElement(i,i);
    }
}

void TH2CB::FillCrystals672(const std::vector<Double_t> &pattern)
{
    if(pattern.size()==672) {

        for(int i=0; i<672; ++i ) {

            SetCrystal672(i, pattern.at(i));
        }
    } else {
        cerr << "TH2CB: FillCrystals672: Wrong pattern size (" << pattern.size() << " / 672)" <<endl;
    }
}

void TH2CB::FillCrystals720(const std::vector<Double_t> &pattern)
{
    if(pattern.size()==720) {

        for(int i=0; i<720; ++i ) {

            SetCrystal720(i,pattern.at(i));

        }
    } else {
        cerr << "TH2CB: FillCrystals720: Wrong pattern size (" << pattern.size() << " / 720)" <<endl;
    }
}

Double_t TH2CB::GetCrystal672(const UInt_t i) const
{
    return GetBinContent(i+1);
}

Double_t TH2CB::GetCrystal720(const UInt_t i) const
{
    const Int_t bin = GetBinOfVBin(i+1);
    if(bin>0)
        return GetBinContent(bin);
    else
        return 0;
}

void TH2CB::SetCrystal672(const UInt_t i, Double_t value)
{
    SetBinContent(i+1, value);
}

void TH2CB::SetCrystal720(const UInt_t i, Double_t value)
{
    const Int_t bin = GetBinOfVBin(i+1);

    if(bin>0)
        return SetBinContent(bin, value);
}

Double_t TH2CB::GetElement(const UInt_t element) const
{
    return GetCrystal720(GetCrystalOfElement(element));
}

void TH2CB::SetElement(const UInt_t element, Double_t value)
{
    SetCrystal720(GetCrystalOfElement(element), value);
}

void TH2CB::SetElements(const std::vector<Double_t> &pattern)
{
    if(pattern.size()==(size_t)GetNumberOfElements()) {

        for(int i=0; i<GetNumberOfElements(); ++i ) {

            SetElement(i,pattern.at(i));

        }
    } else {
        cerr << "TH2CB: FillElements: Wrong pattern size (" << pattern.size() << " / " << GetNumberOfElements() << ")" <<endl;
    }
}

void TH2CB::SetElements(const TH1 &h)
{
    if( h.GetNbinsX() != GetNumberOfElements() ) {
        cerr << "WARNING: Number of bins don't match" << endl;
    }

    const Int_t n = min((Int_t)GetNumberOfElements(), h.GetNbinsX());

    for(Int_t i=1;i<=n; ++i) {
        SetElement(i-1, h.GetBinContent(i));
    }
    SetBinContentChanged(kTRUE);
}

void TH2CB::FillElements(const TH1 &h)
{
    if( h.GetNbinsX() != GetNumberOfElements() ) {
        cerr << "WARNING: Number of bins don't match" << endl;
    }

    const Int_t n = min((Int_t)GetNumberOfElements(), h.GetNbinsX());

    for(Int_t i=1;i<=n; ++i) {
        SetElement(i-1, GetElement(i-1) + h.GetBinContent(i));
    }
    SetBinContentChanged(kTRUE);
}

UInt_t TH2CB::GetCrystalOfElement(const UInt_t element)
{
    if(element<720)
        return crystal_number[element];
    else {
        cerr << "TH2CB::GetCrystalOfElement: Element number out of bounds: (" << element << " / 720)" << endl;
        return 0;
    }
}

UInt_t TH2CB::GetElementOfCrystal(const UInt_t crystal)
{
    if(crystal<720)
        return element_number[crystal];
    else {
        cerr << "TH2CB::GetElementOfCrystal: Crystal number out of bounds: (" << crystal << " / 720)" << endl;
        return 0;
    }
}

void TH2CB::CreateMarker(UInt_t element)
{
    const Int_t bin = GetBinOfVBin(GetCrystalOfElement(element)+1);
    if(bin>0)
        SetMarkerOnBin(bin-1);
}

void TH2CB::SetElements(const TH2CB &h)
{
    TH2Crystals::SetElements(h);
}

void TH2CB::FillElements(const TH2CB &h)
{
    TH2Crystals::FillElements(h);
}

bool TH2CB::IsInHole(const UChar_t a, const UChar_t b, const UChar_t c)
{
    const Int_t vbin = GetVBinOfMMC(a,b,c);
    return bins_in_holes.find(vbin) != bins_in_holes.end();
}

bool TH2CB::IsInHole(const Int_t vbin)
{
    return bins_in_holes.find(vbin) != bins_in_holes.end();
}

/**
 * @brief Puts points for a triable (edge length =1) into a list
 * @return list of points of the triangle
 */
const TH2DrawTool::point_list MakeTriangle() {
    TH2DrawTool::point_list shape(4);
    shape.at(0) = TH2DrawTool::Vector(0.0, 0.0);
    shape.at(1) = TH2DrawTool::Vector(1.0, 0.0);
    shape.at(2) = TH2DrawTool::Vector(0.5, 0.866025);
    shape.at(3) = shape.at(0);

    return shape;
}

const std::set<Int_t> MakeListOfBinsInholes() {

    std::set<Int_t> list;

    // List of Element numbers in holes
    UChar_t holes[2*6*4][3] = {
        {2,1,2},
        {2,1,5},
        {2,1,6},
        {2,1,7},
        {3,2,1},
        {3,2,2},
        {3,2,3},
        {3,2,4},
        {3,3,4},
        {3,3,7},
        {3,3,8},
        {3,3,9},
        {3,1,2},
        {3,1,5},
        {3,1,6},
        {3,1,7},
        {2,2,1},
        {2,2,2},
        {2,2,3},
        {2,2,4},
        {2,3,4},
        {2,3,7},
        {2,3,8},
        {2,3,9},

        {11,1,4},
        {11,1,7},
        {11,1,8},
        {11,1,9},
        {11,3,2},
        {11,3,5},
        {11,3,6},
        {11,3,7},
        {11,4,1},
        {11,4,2},
        {11,4,3},
        {11,4,4},
        {14,1,4},
        {14,1,7},
        {14,1,8},
        {14,1,9},
        {14,3,2},
        {14,3,5},
        {14,3,6},
        {14,3,7},
        {14,4,1},
        {14,4,2},
        {14,4,3},
        {14,4,4}
    };

    // calculate vbin numbers of crystals in holes
    for(UChar_t crystal=0; crystal < 2*6*4; ++crystal) {

        list.insert( TH2CB::GetVBinOfMMC(
                         holes[crystal][0],
                         holes[crystal][1],
                         holes[crystal][2]));

    }

    return list;
}

const std::vector<Int_t> Make_binmap()
{
    std::vector<Int_t> map(721);

    Int_t bin=1;

    for( int major=1;major<=20;++major){
        for( int minor=1;minor<=4;++minor) {
            for( int crystal=1;crystal<=9;++crystal) {

                Int_t vbin = TH2CB::GetVBinOfMMC(major,minor,crystal);


                if( TH2CB::IsInHole(vbin) ) {
                    map[vbin] = -5;
                } else {
                    map[vbin] = bin++;
                }

            }
        }
    }
    return map;
}


