#ifndef TH2TAPS_H
#define TH2TAPS_H

#include "TH2DrawTool.h"
#include "TH2Crystals.h"

namespace ant {

class TH2TAPS: public TH2Crystals {
protected:

    static TH2DrawTool::point_list MakeBaF2Shape();
    static TH2DrawTool::point_list MakePbWO4Shape();

    static const Double_t a;
    static const Double_t b;

    static const TH2DrawTool::point_list baf2_shape;
    static const TH2DrawTool::point_list pbwo4_shape;

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
};

}

#endif
