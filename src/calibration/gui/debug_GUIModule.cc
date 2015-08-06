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

bool DebugModule::Fit(TH1* hist, unsigned channel)
{
    TH2* hist2 = dynamic_cast<TH2*>(hist);

    if(hist2) {
        projection = hist2->ProjectionX("",channel,channel+1);
        func->Fit(projection);

        if(channel==0) {
            return true;
        }

    } else {
        LOG(WARNING) << "Supplied Hist is not 2D";
    }

    // do not show something
    return false;
}

void DebugModule::DisplayFit()
{
    canvas->Show(projection, func.get());
}

void DebugModule::StoreResult(unsigned channel)
{
    LOG(INFO) << "Storing result " << func->GetPeakPosition() << " for channel " << channel;
}

bool DebugModule::Finish()
{
    canvas->Clear();
    canvas->cd();
    TText* text = new TText(0.5, 0.5, "Module finish");
    text->Draw();
    canvas->Draw();
    canvas->Modified();
    canvas->Update();
    return true;
}

void DebugModule::StoreFinish()
{
    LOG(INFO) << "Storing finished result";
}

unsigned DebugModule::GetNumberOfChannels() const
{
    return 10;
}


std::list<CalCanvas*> ant::calibration::gui::DebugModule::GetCanvases() const
{
    return {canvas};
}

void DebugModule::InitGUI()
{
    canvas = new CalCanvas("DebugModule");
}