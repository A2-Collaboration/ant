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
 * @brief Weibull density function
 *
 * f(x)=\lambda \,k\,(\lambda \,x)^{{k-1}}{\mathrm  {e}}^{{-(\lambda \,x)^{k}}}
 * \lambda is the inverse scale parameter and k is the shape parameter
 * Only defined for x >= 0 as well as \lambda > 0 and k > 0
 */
struct weibull {
    static double fct(double* x, double* p);
    static TF1* getTF1();
};

struct landau {
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
 *@brief gaus + polynomial of order n
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
        return helper::makeTF1(Gausexpo::fct, 5 );
    }

};

struct WeibullLandau {
    static double fct(double *x, double *p) {
        return weibull::fct(x, &p[0])*landau::fct(x, &p[2]);
    }

    static TF1* getTF1() {
        return helper::makeTF1(WeibullLandau::fct, 5);
    }
};

/**
 *@brief weibull * landau + polynomial of order n
 */
template <unsigned order>
struct WeibullLandauPol {
    static double fct(double *x, double *p) {
        return weibull::fct(x, &p[0])*landau::fct(x, &p[2]) + pol<order>::fct(x,&p[5]);
    }

    static TF1* getTF1() {
        return helper::makeTF1(WeibullLandauPol<order>::fct, order+6); // weibull*landau: 5 pramams, pol: order + 1
    }
};


}
}
}

