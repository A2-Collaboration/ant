#include "BaseFunctions.h"

#include "base/std_ext.h"

#include <cmath>

#include "TF1.h"

using namespace ant::calibration::functions;

TF1* ant::calibration::functions::helper::makeTF1(double (*fcn)(double *, double *), const unsigned n)
{
    return new TF1("", fcn,0,1,n);
}


double ant::calibration::functions::gaus::fct(double *x, double *p)
{
    return p[0]*exp( - ant::std_ext::sqr(x[0]-p[1]) / (2 * ant::std_ext::sqr(p[2])) );
}

TF1* ant::calibration::functions::gaus::getFT1() {
    return new TF1("",fct,0,1,3);
}


template double pol<2>::fct(double *x, double *p);
template double pol<3>::fct(double *x, double *p);
template double pol<4>::fct(double *x, double *p);

template double GausPol<2>::fct(double *x, double *p);
template double GausPol<3>::fct(double *x, double *p);
template double GausPol<4>::fct(double *x, double *p);


