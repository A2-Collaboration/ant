#include "voigtian.h"

#include <complex>
#include "Faddeeva.hh"

#include "TF1.h"

using namespace ant::math;

using complex = std::complex<double>;

double voigtian::Eval(const double x, const double A, const double mu, const double sigma, const double gamma)
{
    const complex z = (complex(x-mu,0) + complex(0, gamma)) / sigma * M_SQRT1_2;
    const auto w = Faddeeva::w(z);
    const auto v = std::real(w) / sigma * M_SQRT1_2 * 2 * M_2_SQRTPI;
    return A*v;
}

double voigtian::Eval_ROOT(const double* x, const double* p)
{
    return Eval(x[0], p[0], p[1], p[2], p[3]);
}

TF1*voigtian::GetTF1()
{
    TF1* f = new TF1("voigtian", Eval_ROOT, -10, 10, 4);

    f->SetParName(0, "A");
    f->SetParameter(0, 100.0);

    f->SetParName(1, "#mu");
    f->SetParameter(1, 0.0);

    f->SetParName(2, "#sigma");
    f->SetParameter(2, 1.0);

    f->SetParName(3, "#gamma");
    f->SetParameter(3, 1.0);

    return f;
}
