#pragma once

#include "Rtypes.h"
#include "TCanvas.h"

class TFile;
class TTree;

namespace ant {

class TH2CB;
class TH2TAPS;
struct TEvent;

namespace detail {
struct Implementation;
}

struct DisplayClustering : public TCanvas {
protected:
    TH2CB*   h_cb = nullptr;
    TH2TAPS* h_taps = nullptr;
    TTree*   treeEvents = nullptr;
    long long curr_entry = 0;
    TEvent* eventPtr = nullptr;

    void Display();

    // use PIMPL ideom to hide stuff from stupid CINT
    detail::Implementation* impl = nullptr;

public:
    DisplayClustering(TFile* file = nullptr);
    virtual ~DisplayClustering();

    void Next(); //*MENU*
    void Prev(); //*MENU*

    virtual void HandleInput(EEventType button, Int_t x, Int_t y);

    ClassDef(DisplayClustering, 1)
};

}