#pragma once

#include <list>

#ifndef __CINT__
#include "base/interval.h"
#endif

class TF1;


namespace ant {

struct Math {
    static TF1* AsymGaus();

    static TF1* CrystalBall();

};


class TFSum {
protected:
    std::list<TF1*> functions;

    TF1* sum;

    TF1* BuildTF1(const std::string &name, const double min, const double max) const;

public:

    // for ROOT CINT
    TFSum(const std::string &name, TF1* f1, double xmin=-5, double xmax=5);
    TFSum(const std::string &name, TF1* f1, TF1* f2, double xmin=-5, double xmax=5);
    TFSum(const std::string &name, TF1* f1, TF1* f2, TF1* f3, double xmin=-5, double xmax=5);

    //normal ctor
    TFSum(const std::string &name, std::list<TF1*> fs, double xmin=-5, double xmax=5);

    ~TFSum();

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

#ifndef __CINT__
//    interval<double> GetRange() const;
//    void SetRange(const interval<double>& i);
#endif

};

}
