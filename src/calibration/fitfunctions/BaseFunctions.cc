#include "BaseFunctions.h"

#include "base/std_ext/math.h"

#include "TF1.h"
#include "TMath.h"

using namespace std;

using namespace ant::calibration::functions;

TF1* helper::makeTF1(double (*fct)(double *, double *), const unsigned nParameters)
{
    return new TF1("", fct,0,1,nParameters);
}

double gaus::fct(double *x, double *p)
{
    return p[0]*exp( - ant::std_ext::sqr(x[0]-p[1]) / (2 * ant::std_ext::sqr(p[2])) );
}

TF1* gaus::getTF1()
{
    return helper::makeTF1(fct, 3); // new TF1("",fct,0,1,3);
}

double timewalk::fct(double* x, double* p)
{
    if(x[0]<=p[2])
        return 0;

    const auto x0 = x[0]-p[2];

    // using exp and not pow is really needed for numeric stability
    return p[0] + p[5]*x0 +  p[1]*std::exp(-p[4]*x0 - p[3]*std::log(x0));
}

TF1* timewalk::getTF1()
{
    return helper::makeTF1(fct, 6);
}

double exponential::fct(double *x, double *p)
{
    return exp(p[0]+p[1]*x[0]);
}

TF1* exponential::getTF1()
{
    return helper::makeTF1(fct, 3);
}

double landau::fct(double *x, double *p)
{
    return p[0]*TMath::Landau(x[0], p[1], p[2], false);
}

TF1* landau::getTF1()
{
    return helper::makeTF1(fct, 3);
}

double weibull::fct(double *x, double *p)
{
    if (p[0] <= 0 || p[1] <= 0)
        return 0;

    if (x[0] < 0)
        return 0;

    if (x[0] == 0) {
        if (p[1] == 1)
            return p[0];
        else
            return 0;
    }

    if (p[1] == 1)
        return p[0]*exp(-p[0]*x[0]);

    return p[0]*p[1]*std::pow(p[0]*x[0],p[1]-1)*exp(-std::pow(p[0]*x[0],p[1]));
}

TF1* weibull::getTF1()
{
    return helper::makeTF1(fct, 2);
}


template double pol<2>::fct(double *x, double *p);
template double pol<3>::fct(double *x, double *p);
template double pol<4>::fct(double *x, double *p);

template double GausPol<2>::fct(double *x, double *p);
template double GausPol<3>::fct(double *x, double *p);
template double GausPol<4>::fct(double *x, double *p);

template double WeibullLandauPol<1>::fct(double *x, double *p);
template double WeibullLandauPol<2>::fct(double *x, double *p);
