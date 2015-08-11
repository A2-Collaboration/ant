#include "TH2Crystals.h"
#include "TDirectory.h"
#include "TIterator.h"
#include "TGraph.h"
#include "TMarker.h"
#include <iostream>

using namespace ant;
using namespace std;

void TH2Crystals::FillElements(const TH2Crystals &h)
{
    if( h.GetNumberOfBins() == GetNumberOfBins()) {
        TIter thisnext(fBins);
        TIter hnext(h.fBins);

        TH2PolyBin *thisbin=NULL;
        TH2PolyBin *hbin=NULL;

        while(thisnext() && hnext()) {

            thisbin = (TH2PolyBin*) *thisnext;
            hbin = (TH2PolyBin*) *hnext;

            const Double_t v = thisbin->GetContent() + hbin->GetContent();
            thisbin->SetContent(v);

            thisbin->SetChanged(kTRUE);
        }
        SetBinContentChanged(kTRUE);
    } else {
        cerr << "TH2Crystals::FillElements: Number of bins don't match: ( " << GetNumberOfBins() << " / " << h.GetNumberOfBins() << " )" << endl;
    }
}

void TH2Crystals::SetElements(const TH2Crystals &h)
{
    if( h.GetNumberOfBins() == GetNumberOfBins()) {
        TIter thisnext(fBins);
        TIter hnext(h.fBins);

        TH2PolyBin *thisbin=NULL;
        TH2PolyBin *hbin=NULL;

        while(thisnext() && hnext()) {
            thisbin = (TH2PolyBin*) *thisnext;
            hbin = (TH2PolyBin*) *hnext;
            thisbin->SetContent(hbin->GetContent());
            thisbin->SetChanged(kTRUE);
        }
        SetBinContentChanged(kTRUE);
    } else {
        cerr << "TH2Crystals::SetElements: Number of bins don't match: ( " << GetNumberOfBins() << " / " << h.GetNumberOfBins() << " )" << endl;
    }
}

TMarker*TH2Crystals::SetMarkerOnBin(const Int_t bin)
{
    TList* bins = GetBins();
    if(bin < bins->GetEntries()) {
        TH2PolyBin* pbin = dynamic_cast<TH2PolyBin*>(bins->At(bin));
        if(pbin) {
            TGraph* g = dynamic_cast<TGraph*>(pbin->GetPolygon());
            if(g) {
                double x,y;
                calcCOG(g,x,y);
                return new TMarker(x,y,29);
            }
        }
    }
    return nullptr;
}

void TH2Crystals::calcCOG(TGraph* g, double& x, double& y)
{
    if(g->GetN()==0)
        return;

    // save first point
    double first_x,first_y;
    g->GetPoint(0,first_x,first_y);

    x=first_x;
    y=first_y;

    // loop over all points
    double px,py;
    for(Int_t i=1; i<g->GetN();++i) {

        g->GetPoint(i,px,py);

        // and sum up
        x+=px;
        y+=py;
    }

    // if last point == first point (closed polygon): don't sum again/subtract
    if(px==first_x && py == first_y) {
        x -= px;
        y -= py;

        // calc average
        x /= g->GetN()-1;
        y /= g->GetN()-1;
    } else {
        // calc average
        x /= g->GetN();
        y /= g->GetN();
    }
}

TH2Crystals::TH2Crystals( const std::string& name, const std::string& title): TH2Poly()
{
    if(!name.empty())
        SetName(name.c_str());
    if(!title.empty())
        SetTitle(title.c_str());

    SetDirectory(gDirectory);
}

Double_t TH2Crystals::GetElement(const UInt_t element) const
{
    return GetBinContent(element+1);
}

void TH2Crystals::SetElement(const UInt_t element, Double_t value)
{
    SetBinContent( element+1, value);
}

void TH2Crystals::FillBinNumbers() {
    for( Int_t i=0; i < GetNumberOfBins(); ++i ) {
        SetBinContent(i,i);
    }
}

void TH2Crystals::FillElementNumbers()
{
    for( Int_t i=0; i < GetNumberOfElements(); ++i) {
        SetElement(i,i);
    }
}

void TH2Crystals::SetElements(const std::vector<Double_t> &pattern)
{

    if(pattern.size()==(size_t)GetNumberOfElements()) {

        for(size_t i=0; i<pattern.size(); ++i ) {

            SetElement(i,pattern.at(i));

        }
    } else {
        std::cerr << "FillElements: Wrong pattern size (" << pattern.size() << " / " << GetNumberOfElements() << ")" << std::endl;
    }
}

void TH2Crystals::SetElements(const TH1 &h)
{
    if( h.GetNbinsX() != GetNumberOfElements() ) {
        cerr << "WARNING: Number of bis don't match" << endl;
    }

    TIter next(fBins);
    TObject *obj;
    TH2PolyBin *bin;
    Int_t hbin=1;
    while ( (obj = next()) && (hbin <= h.GetNbinsX()) ) {
       bin = (TH2PolyBin*) obj;
       bin->SetContent(h.GetBinContent(hbin++));
    }
}

void TH2Crystals::FillElement(const UInt_t element, const Double_t w)
{
    SetElement( element, GetElement(element) + w);
}

void TH2Crystals::FillElements(const std::vector<Double_t> &pattern)
{
    if(pattern.size()==(size_t)GetNumberOfElements()) {

        for(size_t element=0; element<pattern.size(); ++element ) {

            FillElement(element, pattern.at(element));

        }
    } else {
        std::cerr << "FillElements: Wrong pattern size (" << pattern.size() << " / " << GetNumberOfElements() << ")" << std::endl;
    }
}

void TH2Crystals::FillElements(const TH1 &h)
{
    if( h.GetNbinsX() != GetNumberOfElements() ) {
        cerr << "WARNING: Number of bis don't match" << endl;
    }

    TIter next(fBins);
    TObject *obj;
    TH2PolyBin *bin;
    Int_t hbin=1;
    while ( (obj = next()) && (hbin <= h.GetNbinsX()) ) {
       bin = (TH2PolyBin*) obj;
       bin->SetContent( bin->GetContent() + h.GetBinContent(hbin++));
    }
    SetBinContentChanged(kTRUE);

}

Int_t TH2Crystals::GetNumberOfElements() const
{
    return GetNumberOfBins();
}

void TH2Crystals::ResetElements(const Double_t value)
{
    TIter next(fBins);
    TObject *obj;
    TH2PolyBin *bin;

    while ((obj = next())) {
       bin = (TH2PolyBin*) obj;
       bin->SetContent(value);
    }
    SetBinContentChanged(kTRUE);
}

TMarker*TH2Crystals::SetMarker(const UInt_t element)
{
    return SetMarkerOnBin(element);
}
