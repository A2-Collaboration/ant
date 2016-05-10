#include <fstream>
#include "CB_SourceCalib.h"

#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausexpo.h"

#include "expconfig/detectors/CB.h"
#include "tree/TDetectorReadHit.h"


#include "root-addons/cbtaps_display/TH2CB.h"
#include "base/Logger.h"

#include "TImage.h"


using namespace std;
using namespace ant;
using namespace ant::calibration;

CB_SourceCalib::CB_SourceCalib(std::shared_ptr<expconfig::detector::CB> cb,
                               std::shared_ptr<DataManager> calmgr,
                               Calibration::Converter::ptr_t converter
                               ) :
    Calibration::PhysicsModule("CB_SourceCalib"),
    cb_detector(cb),
    calibrationManager(calmgr),
    Converter(move(converter))

{
}

CB_SourceCalib::~CB_SourceCalib()
{
}

void CB_SourceCalib::ApplyTo(const readhits_t &hits)
{

    const auto& dethits = hits.get_item(Detector_t::Type_t::CB);
    //const auto& dethits = hits.get_item(cb_detector->Type);

    for(TDetectorReadHit& dethit : dethits) {
        if(dethit.ChannelType != Channel_t::Type_t::Integral)
            continue;
        std::vector<double> all_values;
        dethit.Converted = Converter->Convert(dethit.RawData);
        for(double value : dethit.Converted){
            all_values.push_back(value);
        }

        for (double value: all_values)
        {
            dethit.Values.push_back(value);
        }
    }


}
//GUI
void CB_SourceCalib::GetGUIs(list<unique_ptr<gui::CalibModule_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<TheGUI>(
                          GetName(),
                          calibrationManager,
                          cb_detector
                          ));
}

CB_SourceCalib::TheGUI::TheGUI(const string &basename,
                               const std::shared_ptr<DataManager> & calmgr,
                               const std::shared_ptr<expconfig::detector::CB> & cb) :
    gui::CalibModule_traits(basename),
    calibrationManager(calmgr),
    cb_detector(cb),
    func(make_shared<gui::FitGausexpo>())
{
}

CB_SourceCalib::TheGUI::~TheGUI()
{

}

shared_ptr<TH1> CB_SourceCalib::TheGUI::GetHistogram(const WrapTFile &file) const
{
    return file.GetSharedHist<TH1>("CB_SourceCalib/HitsADC_Cluster");
}

unsigned CB_SourceCalib::TheGUI::GetNumberOfChannels() const
{
    return cb_detector->GetNChannels();
}

void CB_SourceCalib::TheGUI::InitGUI(gui::ManagerWindow_traits * window)
{
    canvas = window ->AddCalCanvas();
    h_peaks = new TH1D("h_peaks", "Peak postitions", GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("AmBe Peak");
    AmBe_peaks= new TH1D("AMBe_peaks", "number of Peaks against the ADC channel", 160,0.,160.);

}

void CB_SourceCalib::TheGUI::StartSlice(const interval<TID> & range)
{

}

gui::CalibModule_traits::DoFitReturn_t CB_SourceCalib::TheGUI::DoFit(TH1 *hist, unsigned channel)
{
    if(cb_detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("h_projection",channel+1,channel+1);

    func->SetDefaults(h_projection);
    auto fit_loop = [this] (size_t retries) {
        do {
            func->Fit(h_projection);
            VLOG(8) << "Chi2/dof = " << func->Chi2NDF();
            if(func->Chi2NDF() < AutoStopOnChi2 && func->GetPeakPosition() < 100 && func->GetPeakPosition() >22) {
                return true;
            }
            retries--;
        }
        while(retries>0);
        return false;
    };
    if (fit_loop(8))
        return DoFitReturn_t::Next;

    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF();
    return DoFitReturn_t::Display;

}

void CB_SourceCalib::TheGUI::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void CB_SourceCalib::TheGUI::StoreFit(unsigned channel)
{
    fitParameters[channel] = func->Save();
    const double AmBe_Peak= func->GetPeakPosition();
    const double Peak_Error= func->GetPeakError();
    h_peaks->SetBinContent(channel+1, AmBe_Peak);
    AmBe_peaks->Fill(AmBe_Peak);
    //LOG(INFO) << "Peak-postition:  " << channel  << "  " << AmBe_Peak ;

    fstream file;
    file.open("Peakpositionentest.txt", ios::app);
    file << channel+1 << "\t" << AmBe_Peak << "\t"<< Peak_Error <<endl;
    canvas->Print("Test.pdf(");

//    TImage *img = TImage::Create();
//    img->FromPad(canvas);
//    img->WriteImage("Test22.jpg");

}

bool CB_SourceCalib::TheGUI::FinishSlice()
{
    canvas->Print("Test.pdf)");
    canvas->Clear();
    canvas->Divide(2,1);
    canvas->cd(1);
    h_peaks->SetStats(false);
    h_peaks->Draw();
    canvas->cd(2);
    AmBe_peaks->Draw();


    return true;
}

void CB_SourceCalib::TheGUI::StoreFinishSlice(const interval<TID> &range)
{

}

std::string CB_SourceCalib::TheGUI::GetName() const
{
    return "bla";
}












