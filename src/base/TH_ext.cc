#include "base/TH_ext.h"

#include "TDirectory.h"

namespace ant {

namespace TH_ext {

// copied and adapted from TH2::FitSlicesY/DoFitSlices
TH1D* FitSlicesY(TH2* h, TF1 *f1, Int_t cut, double IQR_range_lo, double IQR_range_hi)
{
    TAxis& outerAxis = *h->GetXaxis();
    Int_t nbins  = outerAxis.GetNbins();

    Int_t npar = f1->GetNpar();

    //Create one histogram for each function parameter

    char *name   = new char[2000];
    char *title  = new char[2000];
    const TArrayD *bins = outerAxis.GetXbins();

    snprintf(name,2000,"%s_Mean",h->GetName());
    snprintf(title,2000,"Fitted value of Mean");
    delete gDirectory->FindObject(name);
    TH1D *hmean = nullptr;
    if (bins->fN == 0) {
        hmean = new TH1D(name,title, nbins, outerAxis.GetXmin(), outerAxis.GetXmax());
    } else {
        hmean = new TH1D(name,title, nbins,bins->fArray);
    }
    hmean->GetXaxis()->SetTitle(outerAxis.GetTitle());

    // Loop on all bins in Y, generate a projection along X
    struct value_t {
        int Bin;
        double Value;
        double Error;
    };
    std::vector<value_t> means;
    for (Int_t bin=0;bin<=nbins+1;bin++) {
        TH1D *hp = h->ProjectionY("_temp",bin,bin,"e");
        if (hp == 0) continue;
        Long64_t nentries = Long64_t(hp->GetEntries());
        if (nentries == 0 || nentries < cut) {delete hp; continue;}


        const double max_pos = hp->GetXaxis()->GetBinCenter(hp->GetMaximumBin());
        const double sigma = hp->GetRMS();

        // setting meaningful start parameters and limits
        // is crucial for a good fit!
        f1->SetParameter(0, hp->GetMaximum());
        f1->SetParLimits(0, 0, hp->GetMaximum());
        f1->SetParLimits(1, max_pos-4*sigma, max_pos+4*sigma);
        f1->SetParameter(1, max_pos);
        f1->SetParLimits(2, 0, 60);
        f1->SetParameter(2, sigma); // set sigma
        f1->SetRange(max_pos-4*sigma, max_pos+4*sigma);

        hp->Fit(f1,"QBNR"); // B important for limits!

        Int_t npfits = f1->GetNumberFitPoints();
        if (npfits > npar && npfits >= cut)
            means.push_back({bin, f1->GetParameter(1), f1->GetParError(1)});

        delete hp;
    }
    delete [] name;
    delete [] title;

    // get some robust estimate of mean errors,
    // to kick out strange outliers
    std_ext::IQR iqr;
    for(const auto& mean : means) {
        iqr.Add(mean.Error);
    }

    auto valid_range = iqr.GetN()<2 ?
                           interval<double>(-std_ext::inf, std_ext::inf) :
                           interval<double>(iqr.GetMedian() - IQR_range_lo*iqr.GetIQR(),
                                            iqr.GetMedian() + IQR_range_hi*iqr.GetIQR());

    for(const auto& mean : means) {
        if(!valid_range.Contains(mean.Error))
            continue;
        hmean->Fill(outerAxis.GetBinCenter(mean.Bin),mean.Value);
        hmean->SetBinError(mean.Bin,mean.Error);
    }
    return hmean;
}

}

}
