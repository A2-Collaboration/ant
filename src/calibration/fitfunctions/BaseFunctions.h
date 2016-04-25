#pragma once

class TF1;

namespace ant {
namespace calibration {
namespace functions {

struct helper {
    static TF1* makeTF1(double (*fcn)(double*, double*), const unsigned nParameters);
};

struct gaus {
    static double fct(double* x, double* p);
    static TF1* getTF1();
};

struct timewalk {
    static double fct(double* x, double* p);
    static TF1* getTF1();
};

struct exponential {
    static double fct(double* x, double* p);
    static TF1* getTF1();
};


/**
 *@brief Polynomial of order n
 */
template <unsigned order>
struct pol {
    static double fct(double* x, double* p) {

        double res  = p[0];
        double mult = 1;

        for(unsigned i=1; i<=order; ++i) {
            mult *= x[0];
            res += p[i] * mult;
        }

        return res;
    }

    static TF1* getTF1() {
        return helper::makeTF1(pol<order>::fct,order+1);
    }
};

/**
 *@brief gaus + polynomila of order n
 */
template <unsigned order>
struct GausPol {
    static double fct(double *x, double *p) {
        return gaus::fct(x, &p[0]) + pol<order>::fct(x,&p[3]);
    }

    static TF1* getTF1() {
        return helper::makeTF1(GausPol<order>::fct,order+4); // gaus: 3 pramams, pol: order + 1
    }
};

struct Gausexpo {
    static double fct(double *x, double *p) {
        return gaus::fct(x, &p[0])+exponential::fct(x, &p[3]);
    }

    static TF1* getTF1() {
        return helper::makeTF1(Gausexpo::fct, 6 );
    }

};


}
}
}

