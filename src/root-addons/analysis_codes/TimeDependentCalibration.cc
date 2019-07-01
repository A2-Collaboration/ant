#include "TimeDependentCalibration.h"

#include "base/WrapTFile.h"
#include "analysis/plot/HistogramFactory.h"
#include "tree/TAntHeader.h"
#include "expconfig/ExpConfig.h"

#include "TH2D.h"
#include "TF1.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;

void TimeDependentCalibration::MakeCBEnergyFile(const char* basefilename,
                                                const char* setupname,
                                                int fillsPerChannel, int nSlices)
{
    vector<double> peakPos(nSlices, 135);

    if(nSlices>20) {
        int i = 10;
        peakPos[i++] = 138;
        peakPos[i++] = 140;
        peakPos[i++] = 137;
        peakPos[i++] = 133;
    }

    if(nSlices>50) {
        int i = 30;
        peakPos[i++] = 136;
        peakPos[i++] = 137;
        peakPos[i++] = 138;
        peakPos[i++] = 139;
        peakPos[i++] = 140;
        peakPos[i++] = 141;
        peakPos[i++] = 139;
        peakPos[i++] = 137;
        peakPos[i++] = 135;
        peakPos[i++] = 133;
        peakPos[i++] = 131;
        peakPos[i++] = 129;
    }

    ExpConfig::Setup::SetByName(setupname);
    auto cb = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);
    const auto nCBChannels = cb->GetNChannels();

    for(auto slice=0;slice<nSlices;slice++) {

        WrapTFileOutput out(std_ext::formatter() << basefilename << "_" << slice << ".root", true);
        TAntHeader* header = new TAntHeader();
        gDirectory->Add(header);

        HistogramFactory HistFac("CB_Energy");
        auto ggIM = HistFac.makeTH2D("ggIM",{"IM / MeV",1000},{"Channel",nCBChannels},"ggIM");

        TF1 f("f", "gaus+expo(3)*pol2(5)");
        f.SetRange(0, 1000);
        f.SetNpx(1000);
        f.SetParameter(0, 1000);
        f.SetParameter(1,  peakPos.at(slice));
        f.SetParameter(2,    8);
        f.SetParameter(3,  3); // expo = exp(p0+p1*x)
        f.SetParameter(4,  -0.01);
        f.SetParameter(5,  0);
        f.SetParameter(6,  1);
        f.SetParameter(7,  0);

        for(auto ch=0u;ch<nCBChannels;ch++) {
            for(auto i=0;i<fillsPerChannel;i++) {
                ggIM->Fill(f.GetRandom(), ch);
            }
        }

        header->CmdLine = "root ant::TimeDependentCalibration";
        header->FirstID = TID(slice, 0, {TID::Flags_t::AdHoc});
        header->LastID = TID(slice, 1, {TID::Flags_t::AdHoc});
        header->SetupName = setupname;

    }
}
