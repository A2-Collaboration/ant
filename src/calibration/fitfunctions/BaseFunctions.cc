#include "BaseFunctions.h"

#include "base/std_ext.h"

#include <cmath>

#include "TF1.h"

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
    return p[0] + p[1]/std::pow(x[0] + p[2], p[3]);
}

TF1* timewalk::getTF1()
{
    return helper::makeTF1(fct, 4);
}


template double pol<2>::fct(double *x, double *p);
template double pol<3>::fct(double *x, double *p);
template double pol<4>::fct(double *x, double *p);

template double GausPol<2>::fct(double *x, double *p);
template double GausPol<3>::fct(double *x, double *p);
template double GausPol<4>::fct(double *x, double *p);
