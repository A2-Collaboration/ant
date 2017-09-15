#include "SigmaPlus_tools.h"

#include <numeric>

using namespace ant;
using namespace std;
using namespace RooFit;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;


namespace SIGMA {

LumiFitter_t::result_t LumiFitter_t::DoFit(TH1D* histLumi)
{
    result_t result;
    
    for (auto i = 1 ; i < histLumi->GetNbinsX() ; ++i)
    {
        result.Fitted_Egs.emplace_back(histLumi->GetBinCenter(i),histLumi->GetBinWidth(i) / sqrt(12.));
        result.Fitted_Lumis.emplace_back(histLumi->GetBinContent(i),histLumi->GetBinError(i));
//        cout << result.Fitted_Egs.back() << endl;
//        cout << result.Fitted_Lumis.back() << endl << endl;
    }
    
    auto residuals = [] (const Value_t& I,
                     const Value_t& c1, const Value_t& c2, const Value_t& c3, const Value_t& c4, const Value_t& c5,
                     const vector<Value_t>& Es, const vector<Value_t>& Ls)
    {
        vector<double> residuals(Ls.size());
        transform(Es.begin(),Es.end(),Ls.begin(),residuals.begin(),
                  [&I,&c1,&c2,&c3,&c4,&c5](const double& Es_i, const double& Ls_i)
        {

            return I
                    + c1 * Es_i
                    + c2 * pow(Es_i,2)
                    + c3 * pow(Es_i,3)
                    + c4 * pow(Es_i,4)
//                    + c5 * pow(Es_i,5)   WTF??????!!!!!!
                    - Ls_i;
        });
        return  residuals;
    };
    
    fitter.DoFit(result.Intercept, result.c1,result.c2,result.c3,result.c4,result.c5,
                 result.Fitted_Egs,result.Fitted_Lumis,residuals);
    
    return result;
}


}

