#include "debug_GUIModule.h"
#include "GUIbase.h"
#include "CalCanvas.h"
#include "FitGausPol3.h"
#include "base/Logger.h"
#include "calibration/gui/FitGausPol3.h"

#include "TH2.h"
#include "TText.h"

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;


DebugModule::DebugModule() :
    func(make_shared<FitGausPol3>())
{
}

DebugModule::~DebugModule()
{

}

string DebugModule::GetHistogramName() const
{
    return "Calibration_CB_Energy_Gains/ggIM";
}

GUIClientInterface::FitStatus DebugModule::Fit(CalCanvas* c, TH1* hist, unsigned channel)
{
    TH2* hist2 = dynamic_cast<TH2*>(hist);

    if(hist2) {
        TH1* hist1 = hist2->ProjectionX("",channel,channel+1);
        func->Fit(hist1);

        if(channel==0) {
            c->Show(hist1, func.get());
            return FitStatus::GUIWait;
        }

    } else {
        LOG(WARNING) << "Supplied Hist is not 2D";
    }

    return FitStatus::FitOK;
}

void DebugModule::StoreResult(unsigned channel)
{
    LOG(INFO) << "Storing result " << func->GetPeakPosition() << " for channel " << channel;
}

GUIClientInterface::FitStatus DebugModule::Finish(CalCanvas* c)
{
    c->cd();
    TText* text = new TText(0.5, 0.5, "Module finish");
    text->Draw();
    c->Draw();
    c->Modified();
    c->Update();
    return FitStatus::GUIWait;
}

void DebugModule::StoreFinish()
{
    LOG(INFO) << "Storing finished result";
}

unsigned DebugModule::GetNumberOfChannels()
{
    return 720;
}
