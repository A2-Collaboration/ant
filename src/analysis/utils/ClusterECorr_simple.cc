#include "ClusterECorr_simple.h"


using namespace std;
using namespace ant;
using namespace ant::analysis::utils;
using namespace ant::std_ext;

ClusterECorr_simple::ClusterECorr_simple() {}

ClusterECorr_simple::~ClusterECorr_simple() {}

void ClusterECorr_simple::LoadECorr(const string& filename, const string& hname)
{
    TFile f(Form("%s",filename.c_str()));
    if(f.IsZombie()){
        return;
    }
    else {
        TH1D *hOrig = dynamic_cast<TH1D*>(f.Get(Form("%s",hname.c_str())));
        if(!hOrig)
            throw(std::runtime_error("Histogram not found: "+hname));
        else
            hECorrSet = true;
        hECorr = *dynamic_cast<TH1D*>(hOrig->Clone());
        hOrig->Delete();
    }
}

double ClusterECorr_simple::GetECorr(const double CluEin)
{
    if(hECorrSet){
        double CluE = CluEin;
        double CluEMin = hECorr.GetXaxis()->GetBinCenter(1);
        int maxbin = hECorr.GetXaxis()->GetNbins();
        double CluEMax = hECorr.GetXaxis()->GetBinCenter(maxbin);
        if(CluE<CluEMin) CluE = CluEMin;
        if(CluE>CluEMax) CluE = CluEMax;
        int bin = hECorr.FindFixBin(CluE);
        double ECorr = hECorr.GetBinContent(bin);
        return ECorr;
    }
    else
        throw(std::runtime_error("Correction histogram not set"));
    return 0;
}
