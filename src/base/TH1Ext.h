#pragma once

#include "TH1.h"

namespace ant{

double GetMaxPos(TH1* hist){
    return hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
}

}
