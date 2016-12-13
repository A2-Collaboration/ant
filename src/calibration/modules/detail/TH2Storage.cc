#include "TH2Storage.h"
#include "tree/TCalibrationData.h"
#include "TH2.h"
#include "TH2D.h"

using namespace ant;
using namespace ant::calibration::detail;

void TH2Storage::Encode(const TH2* histogram, ant::TCalibrationData& cdata)
{
    cdata.FitParameters.clear();
    cdata.Data.clear();

    cdata.FitParameters.emplace_back(TCalibrationData::TFitParameters(0,{double(histogram->GetNbinsX()), histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax()}));
    cdata.FitParameters.emplace_back(TCalibrationData::TFitParameters(1,{double(histogram->GetNbinsY()), histogram->GetYaxis()->GetXmin(), histogram->GetYaxis()->GetXmax()}));

   cdata.Data.reserve(size_t(histogram->GetNbinsX()*histogram->GetNbinsY()));

   for(int y=1; y<=histogram->GetNbinsY(); ++y) {
       for(int x=1; x <= histogram->GetNbinsX(); ++x) {
           cdata.Data.emplace_back(TCalibrationData::Entry(0, histogram->GetBinContent(x,y)));
       }
   }
}

TH2D* TH2Storage::Decode(const ant::TCalibrationData& cdata)
{
    const auto xb = GetXBins(cdata);
    const auto yb = GetYBins(cdata);

    auto hist = new TH2D("",cdata.CalibrationID.c_str(),
                         int(xb.Bins()), xb.Start(), xb.Stop(),
                         int(yb.Bins()), yb.Start(), yb.Stop());

    auto p = cdata.Data.cbegin();
    for(int y=1; y<=int(yb.Bins()); ++y) {
        for(int x=1; x <= int(xb.Bins()); ++x) {
            hist->SetBinContent(x,y,p->Value);
            ++p;
        }
    }

    return hist;
}

BinSettings TH2Storage::GetXBins(const TCalibrationData& cdata)
{
    const auto& xp = cdata.FitParameters.at(0).Value;
    return {unsigned(xp.at(0)), xp.at(1), xp.at(2)};
}

BinSettings TH2Storage::GetYBins(const TCalibrationData& cdata)
{
    const auto& yp = cdata.FitParameters.at(1).Value;
    return {unsigned(yp.at(0)), yp.at(1), yp.at(2)};
}
