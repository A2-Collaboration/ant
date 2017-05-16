#pragma once

#include <list>

#ifndef __CINT__
#include "base/interval.h"
#include "base/piecewise_interval.h"
#include "analysis/plot/RootDraw.h"
#endif

#include "TF1.h"

class TH1;
class TFitResultPtr;

namespace ant {

struct Math {
    static TF1* AsymGaus();

    static TF1* CrystalBall();

};

#ifndef __CINT__
class TFSum : public root_drawable_traits
#else
class TFSum
#endif
{

protected:
    std::list<TF1*> functions;

    TF1* sum;

    TF1* BuildTF1(const std::string &name, const double min, const double max) const;

    static void syncTF1Par(TF1* src, TF1* dst);
    static void syncTF1Par(TF1* src, TF1* dst, int src_j, int dst_j);
public:

    // for ROOT CINT
    TFSum(const std::string &name, TF1* f1, double xmin=-5, double xmax=5);
    TFSum(const std::string &name, TF1* f1, TF1* f2, double xmin=-5, double xmax=5);
    TFSum(const std::string &name, TF1* f1, TF1* f2, TF1* f3, double xmin=-5, double xmax=5);

    //normal ctor
    TFSum(const std::string &name, std::list<TF1*> fs, double xmin=-5, double xmax=5);

    virtual ~TFSum();

    /**
     * @brief Evalute the function
     * @param x
     * @param p
     * @return
     */
    double operator() (const double* x, const double* p) const;

    unsigned GetNpar() const;

    /**
     * @brief Copy Parameters, Errors and Limits to the Sum TF1
     */
    void SyncToSum();

    /**
     * @brief Copy Parameters, Errors, and Limits to the individual functions
     */
    void SyncToFcts();

    void Draw();

    void SetRange(double xmin, double xmax);

    void GetRange(double& xmin, double &xmax) const;

    TF1* Function() { return sum; }

    operator TF1*() {
        return sum;
    }

    void SetNpx(int n);


    static TFitResultPtr FitRanged(TH1* h, TF1* f,
                                   double x_low, double x_high,
                                   const std::string& fitopts = "REM0NB");
    static TFitResultPtr FitRanged(TH1* h, TF1* f,
                                   double x1_low, double x1_high, double x2_low, double x2_high,
                                   const std::string& fitopts = "REM0NB");

    static TF1* MakeRanged(TF1* f, double x_low, double x_high);
    static TF1* MakeRanged(TF1* f, double x1_low, double x1_high, double x2_low, double x2_high);




#ifndef __CINT__
    static TFitResultPtr FitRanged(
            TH1* h, TF1* f,
            const PiecewiseInterval<double>& range,
            const std::string& fitopts = "REM0NB");

    static TF1* MakeRanged(TF1* f, const PiecewiseInterval<double>& range);

    template<typename Filter>
    static TF1* MakeFiltered(TF1* f, Filter filter) {
        auto filtered_fct = [f, filter] (double* x, double* p) {
            if(!filter(x))
                TF1::RejectPoint();
            return f->EvalPar(x, p);
        };
        auto f_filtered = new TF1((f->GetName()+std::string("_filtered")).c_str(),
                                  filtered_fct,
                                  f->GetXmin(), f->GetXmax(), f->GetNpar());
        syncTF1Par(f, f_filtered);
        return f_filtered;
    }

    virtual void Draw(const std::string& option) const;
#endif

};

}
