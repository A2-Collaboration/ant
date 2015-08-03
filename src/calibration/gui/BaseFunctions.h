#pragma once

class TF1;

namespace ant {
namespace calibration {
namespace functions {

struct helper {
    static TF1* makeTF1(double (*fcn)(double *, double *), const unsigned n);
};

struct gaus {
    static double fct(double *x, double *p);
    static TF1* getFT1();
};


/**
 *@brief Polynomial of order n
 */
template <unsigned order>
struct pol {
    static double fct(double* x, double* p) {

        double res  = 0;
        double mult = 1;

        for(unsigned i=0; i<order; ++i) {
            res += p[i] * mult;
            mult *= x[0];
        }

        return res;
    }

    static TF1* getTF1() {
        return helper::makeTF1(pol<order>::fct,order);
    }
};

/**
 *@brief gaus + polynomila of order n
 */
template <unsigned order>
struct GausPol {
    static double fct(double *x, double *p) {
        return ant::calibration::functions::gaus::fct(x, &p[0]) + pol<order>::fct(x,&p[3]);
    }

    static TF1* getTF1() {
        return helper::makeTF1(GausPol<order>::fct,order+3);
    }
};

}
}
}

