#include "base/TH_ext.h"

#include "TDirectory.h"
#include "base/std_ext/string.h"
#include <iomanip>

using namespace std;

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

TH1D *HistOfBins(const TH2 *h2)
{
    const auto nx = h2->GetNbinsX();
    const auto ny = h2->GetNbinsY();

    auto h = new TH1D("","",  int(sqrt(nx*ny)), h2->GetMinimum(), h2->GetMaximum() + std::numeric_limits<double>::epsilon());
    h->SetXTitle(h2->GetZaxis()->GetTitle());
    for(int x=1; x<=nx; ++x) {
        for(int y=1; y<=ny; ++y) {
            h->Fill(h2->GetBinContent(x,y));
        }
    }
    return h;
}

std::string TH1ToLaTeX(const TH1 *h, const int precission)
{
    std_ext::formatter f;

    f << "\\begin{tabular}{rr}\n"
      << "\t\\textbf{" << h->GetXaxis()->GetTitle() << "} & \\textbf{" << h->GetYaxis()->GetTitle() << "} \\\\\n";
    for(int i=1;i<=h->GetNbinsX();++i) {
        f << "\t\\num{" << setw(8) << setprecision(precission) << h->GetBinCenter(i) << "} & "
          << "\\num{" << setw(8) << setprecision(precission) << h->GetBinContent(i) << "} \\\\\n";
    }
    f << "\\end{tabular}";

    return f;
}

TH2D* GetSlice(const TH3& h, const int b, const char* projection)
{
    string project = std_ext::to_lower(projection);
    const TAxis *sliceAxis, *axis1, *axis2;

    if (std_ext::contains(project, "xy")) {
        sliceAxis = h.GetZaxis();
        axis1 = h.GetYaxis();
        axis2 = h.GetXaxis();
    } else if (std_ext::contains(project, "xz")) {
        sliceAxis = h.GetYaxis();
        axis1 = h.GetZaxis();
        axis2 = h.GetXaxis();
    } else if (std_ext::contains(project, "yx")) {
        sliceAxis = h.GetZaxis();
        axis1 = h.GetXaxis();
        axis2 = h.GetYaxis();
    } else if (std_ext::contains(project, "yz")) {
        sliceAxis = h.GetXaxis();
        axis1 = h.GetZaxis();
        axis2 = h.GetYaxis();
    } else if (std_ext::contains(project, "zx")) {
        sliceAxis = h.GetYaxis();
        axis1 = h.GetXaxis();
        axis2 = h.GetZaxis();
    } else if (std_ext::contains(project, "zy")) {
        sliceAxis = h.GetXaxis();
        axis1 = h.GetYaxis();
        axis2 = h.GetZaxis();
    } else {
        cerr << "ERROR: Invalid projection choice! Argument has to be "
             << "a combination of axes, like \"yx\" or \"xz\"" << endl;
        return nullptr;
    }

    int bin = b;
    int n;

    if (bin < 1) {
        cerr << "WARNING: Chosen bin is smaller than axis bins, use bin 1" << endl;
        bin = 1;
    } else if (bin > (n = sliceAxis->GetNbins())) {
        cerr << "WARNING: Chosen bin is greater than number of bins, use bin " << n << endl;
        bin = n;
    }

    stringstream name, title;
    name << h.GetName() << "_project_" << project << "_bin" << bin;
    title << h.GetTitle() << " " << project << " projection bin " << bin;

    return GetSlice(h, name.str().c_str(), title.str().c_str(), bin, axis1, axis2);
}

TH2D* GetSlice(const TH3& hist, const char* name, const char* title, const int bin,
               const TAxis* const axis1, const TAxis* const axis2)
{
    TH2D* res = nullptr;

    if (!axis1 || !axis2) {
        cerr << "Need two axes of the 3D histogram to obtain a slice" << endl;
        return res;
    }

    if (axis1->GetNbins() <= 0 || axis2->GetNbins() <= 0) {
        cerr << "Histogram axes need to have bins" << endl;
        return res;
    }

    TObject* o = gDirectory->FindObject(name);
    if (o)
        delete o;

    // link reference pointers to bin counters to match the correct axes while retrieving bin contents
    int binx, biny, binz = bin;
    int *refx, *refy, *refz;

    if (axis1 == hist.GetXaxis() && axis2 == hist.GetYaxis()) { refx = &binx; refy = &biny; refz = &binz; }  // yx
    if (axis1 == hist.GetYaxis() && axis2 == hist.GetXaxis()) { refx = &biny; refy = &binx; refz = &binz; }  // xy
    if (axis1 == hist.GetXaxis() && axis2 == hist.GetZaxis()) { refx = &binx; refy = &binz; refz = &biny; }  // zx
    if (axis1 == hist.GetZaxis() && axis2 == hist.GetXaxis()) { refx = &biny; refy = &binz; refz = &binx; }  // xz
    if (axis1 == hist.GetYaxis() && axis2 == hist.GetZaxis()) { refx = &binz; refy = &binx; refz = &biny; }  // zy
    if (axis1 == hist.GetZaxis() && axis2 == hist.GetYaxis()) { refx = &binz; refy = &biny; refz = &binx; }  // yz

    res = new TH2D(name, title,
                   axis1->GetNbins(), axis1->GetBinLowEdge(axis1->GetFirst()), axis1->GetBinUpEdge(axis1->GetLast()),
                   axis2->GetNbins(), axis2->GetBinLowEdge(axis2->GetFirst()), axis2->GetBinUpEdge(axis2->GetLast()));
    res->GetXaxis()->ImportAttributes(axis1);
    res->GetYaxis()->ImportAttributes(axis2);

    for (binx = 0; binx <= res->GetNbinsX()+1; ++binx)
        for (biny = 0; biny <= res->GetNbinsY()+1; ++biny)
            res->SetBinContent(binx, biny, hist.GetBinContent(*refx,*refy,*refz));

    res->ResetStats();

    return res;
}

}

}
