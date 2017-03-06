#include "TimeDependentCalibration.h"

#include "base/WrapTFile.h"
#include "analysis/plot/HistogramFactory.h"

#include "TH2D.h"
#include "TF1.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;

void TimeDependentCalibration::MakeCBEnergyFile(const char* outfilename, int fillsPerChannel)
{
    constexpr auto nCBChannels = 720; // too lazy to use ExpConfig here :)
    WrapTFileOutput out(outfilename, WrapTFileOutput::mode_t::recreate, true);
    HistogramFactory HistFac("CB_Energy");
    auto ggIM = HistFac.makeTH2D("ggIM",{"IM / MeV",{1000}},{"Channel",{nCBChannels}},"ggIM");
//    ggIM->Draw();

    auto f = new TF1("f", "gaus+expo(3)*pol2(5)");
    f->SetRange(0, 1000);
    f->SetNpx(1000);
    f->SetParameter(0, 1000);
    f->SetParameter(1,  135);
    f->SetParameter(2,    8);
    f->SetParameter(3,  3); // expo = exp(p0+p1*x)
    f->SetParameter(4,  -0.01);
    f->SetParameter(5,  0);
    f->SetParameter(6,  1);
    f->SetParameter(7,  0);

    for(auto ch=0;ch<nCBChannels;ch++) {
        for(auto i=0;i<fillsPerChannel;i++) {
            ggIM->Fill(f->GetRandom(), ch);
        }
    }

    //    ggIM->Draw("colz");
//    f->Draw();
}
