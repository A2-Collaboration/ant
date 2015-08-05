#include "debug_GUIModule.h"
#include "GUIbase.h"
#include "FitCanvas.h"
#include "FitGausPol3.h"
#include "base/Logger.h"
#include "calibration/gui/FitGausPol3.h"

#include "TH2.h"

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;


DebugModule::DebugModule(): func(new FitGausPol3())
{
}

DebugModule::~DebugModule()
{

}

string DebugModule::GetHistogramName() const
{
    return "Calibration_CB_Energy_Gains/ggIM";
}

GUIClientInterface::FitStatus DebugModule::Fit(TH1* hist, unsigned channel)
{
    TH2* hist2 = dynamic_cast<TH2*>(hist);

    if(hist2) {
        TH1* hist1 = hist2->ProjectionX("",channel,channel+1);
        func->Fit(hist1);

        if(channel==0) {
            CalCanvas* c = new CalCanvas("aaa");
            c->Draw();
            c->Show(hist1,func);
            guiinstance = c;
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
    guiinstance = nullptr;
}

GUIClientInterface::FitStatus DebugModule::Finish()
{
    return FitStatus::FitOK;
}

TQObject* DebugModule::GetGUIInstance()
{
    return guiinstance;
}

unsigned DebugModule::GetNumberOfChannels()
{
    return 720;
}
