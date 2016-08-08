#ifndef BROWSEHISTOGRAMS_H
#define BROWSEHISTOGRAMS_H

#include <list>

class TH1;

#include "TDirectory.h"
#include "TCanvas.h"

namespace ant {

class BrowseHistogramsCanvas : public TCanvas {
protected:
    typedef std::list<TH1*> hlist;
    hlist hists;
    hlist::const_iterator current;

    std::string draw_option;

public:
    BrowseHistogramsCanvas(TDirectory* dir=gDirectory);
    virtual ~BrowseHistogramsCanvas();

    virtual bool Next(); //*MENU*
    virtual bool Prev(); //*MENU*
    virtual void Scan(TDirectory* dir=gDirectory);  //*MENU*
    virtual void SetOption(const std::string& option); //*MENU*

    virtual void DrawCurrent();

    virtual void HandleInput(EEventType button, Int_t x, Int_t y);

    ClassDef(BrowseHistogramsCanvas, 1)
};

}

#endif
