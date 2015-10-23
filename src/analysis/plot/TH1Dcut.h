#pragma once

#include <TH1D.h>

namespace ant {
namespace analysis {

class TH1Dcut : public TH1D {
    using TH1D::TH1D;

    virtual Int_t Fill(Double_t x) override;

};

}
}


