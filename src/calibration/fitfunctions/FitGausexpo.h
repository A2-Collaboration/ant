#pragma once

#include "FitFunction.h"
#include "TSpectrum.h"



namespace ant{
namespace calibration{
namespace gui {

class FitGausexpo: public PeakingFitFunction {
protected:
    TF1* signal;
    TF1* bg;



public:
    FitGausexpo();

    virtual ~FitGausexpo();

    void Draw() override;
    void Fit(TH1* hist) override;
    void FitSignal(TH1 * hist) override;
    void FitBackground(TH1 * hist) override;
    virtual void SetDefaults(TH1* hist) override;
    void SetRange(ant::interval<double> i) override;
    ant::interval<double> GetRange() const override;
    virtual void Sync() override;
    std::vector<double> Save() const override;
    void Load(const std::vector<double> &data) override;

    virtual double GetPeakPosition() const override;
    virtual double GetPeakWidth() const override;





};

}
}
}
