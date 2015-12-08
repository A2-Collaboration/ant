#ifndef TH2TAPS_H
#define TH2TAPS_H

#include "TH2Crystals.h"

namespace ant {

class TH2DrawTool;

class TH2TAPS: public TH2Crystals {
protected:

    void DrawShape(TH2DrawTool &c, bool isBaF2 = false);

    virtual void Build();

public:
    TH2TAPS( const std::string& name="", const std::string& title="");
    virtual ~TH2TAPS() {}

    // pulling in methods from base class
    using TH2Crystals::SetElements;
    using TH2Crystals::FillElements;
    virtual void FillElements( const TH2TAPS& h);
    virtual void SetElements( const TH2TAPS& h);

    ClassDef(TH2TAPS,2);
};

}

#endif
