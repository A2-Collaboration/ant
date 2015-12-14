#include "TAPSEnergy_Check.h"


#include "expconfig/ExpConfig.h"
#include "analysis/plot/root_draw.h"
#include "base/interval.h"
#include "base/std_ext/math.h"

#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"


#include <iostream>
#include <sstream>

using namespace std;
using namespace ant;



void TAPSEnergy_Check::Analyse(TFile* file)
{
    auto h = dynamic_cast<TH3D*>(file->Get("TAPS_Energy/ggIM_mult"));

    canvas c;
    for(int m=0;m<=10;m++) {
        h->GetZaxis()->SetRangeUser(m,m+1);
        stringstream ss_name;
        ss_name << "Mult_" << m << "_yx";
        c << drawoption("colz") << h->Project3D(ss_name.str().c_str());

    }
    c << dynamic_cast<TH2D*>(file->Get("TAPS_Energy/ggIM"));
    c << endc;
}

