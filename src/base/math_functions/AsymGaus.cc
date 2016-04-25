#include "AsymGaus.h"

#include "TF1.h"
#include "TMath.h"

using namespace ant::math;


double AsymGaus::Eval(const double x, const double A, const double mu, const double sigma_l, const double sigma_r)
{
    return A*TMath::Gaus(x,mu, x>mu ? sigma_r : sigma_l, false);
}

double AsymGaus::Eval_ROOT(const double* x, const double* p)
{
    return Eval(x[0], p[0], p[1], p[2], p[3]);
}

TF1* AsymGaus::GetTF1()
{
    TF1* f = new TF1("asymgaus", Eval_ROOT, -10, 10, 4);
    f->SetNpx(1000);

    f->SetParName(0, "A");
    f->SetParameter(0, 100.0);

    f->SetParName(1, "#mu");
    f->SetParameter(1, 0.0);

    f->SetParName(2, "#sigma_{left}");
    f->SetParameter(2, 1.0);

    f->SetParName(3, "#sigma_{right}");
    f->SetParameter(3, 1.0);

    return f;
}
