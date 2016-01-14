#include "BrowseHistograms.h"
#include "TH1.h"
#include "TKey.h"
#include "TDirectory.h"

#include <iostream>

using namespace std;


BrowseHistogramsCanvas::BrowseHistogramsCanvas(TDirectory* dir):
    draw_option("colz")
{
    Scan(dir);
}

BrowseHistogramsCanvas::~BrowseHistogramsCanvas()
{}

bool BrowseHistogramsCanvas::Next()
{
    if(++current != hists.end()) {
        DrawCurrent();
        return true;
    }
    return false;
}

bool BrowseHistogramsCanvas::Prev()
{
    if(current != hists.begin()) {
        --current;
        DrawCurrent();
        return true;
    }
    return false;
}

void BrowseHistogramsCanvas::Scan(TDirectory* dir)
{

        TList* keys = dir->GetListOfKeys();

        if(keys) {

        TH1* objectPtr = nullptr;
        TKey* key = nullptr;
        TIter nextk(keys);

        while( (key = dynamic_cast<TKey*>(nextk())) )
        {
            objectPtr = dynamic_cast<TH1*>( key->ReadObj() );

            if ( !objectPtr )
                continue;

            hists.push_back(objectPtr);
        }
    }

    current = hists.cbegin();

    DrawCurrent();
}

void BrowseHistogramsCanvas::SetOption(const std::string& option)
{
    draw_option = option;
}

void BrowseHistogramsCanvas::DrawCurrent()
{
    this->cd();
    (*current)->Draw(draw_option.c_str());
    Modified();
    Update();
}

void BrowseHistogramsCanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    if(button == kKeyPress) {
        if(char(x) == 'd') {
            Next();
        } else if( char(x) == 'a') {
            Prev();
        }
    }
    TCanvas::HandleInput(button, x ,y);
}
