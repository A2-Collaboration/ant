#include "Math.h"
#include "base/math_functions/AsymGaus.h"
#include "base/math_functions/CrystalBall.h"

#include "TF1.h"


TF1*ant::Math::AsymGaus()
{
    return math::AsymGaus::GetTF1();
}

TF1 *ant::Math::CrystalBall()
{
    return math::CrystalBall::GetTF1();
}

void ant::TFSum::SyncToSum()
{
    unsigned par = 0;
    for(const auto& f : functions) {
        for(int j=0; j<f->GetNpar(); ++j) {
            sum->SetParameter(par, f->GetParameter(j));
            double low, high;
            f->GetParLimits(j, low, high);
            sum->SetParLimits(par, low, high);
            sum->SetParError(par, f->GetParError(j));
            sum->SetParName(par, f->GetParName(j));
            ++par;
        }
    }
}

void ant::TFSum::SyncToFcts()
{
    unsigned par = 0;
    for(auto& f : functions) {
        for(int j=0; j<f->GetNpar(); ++j) {
            f->SetParameter(j, sum->GetParameter(par));
            double low, high;
            sum->GetParLimits(par, low, high);
            f->SetParLimits(j, low, high);
            f->SetParError(j, sum->GetParError(par));
            f->SetParName(j, sum->GetParName(par));
            ++par;
        }
    }
}

void ant::TFSum::Draw()
{
    sum->Draw("same");
    for(auto& f : functions) {
        f->Draw("same");
    }
}

TF1 *ant::TFSum::BuildTF1(const std::string& name, const double min, const double max) const
{
    unsigned pars = 0;
    for(const auto& f : functions) {
        pars += f->GetNpar();
    }

    return new TF1(name.c_str(), this, min, max, pars);
}

ant::TFSum::TFSum(const std::string &name, TF1 *f1, double xmin, double xmax):
    TFSum(name, std::list<TF1*>({f1}), xmin, xmax)
{}

ant::TFSum::TFSum(const std::string &name, TF1 *f1, TF1 *f2, double xmin, double xmax):
    TFSum(name, {f1,f2}, xmin, xmax)
{}

ant::TFSum::TFSum(const std::string &name, TF1 *f1, TF1 *f2, TF1* f3, double xmin, double xmax):
    TFSum(name, {f1,f2,f3}, xmin, xmax)
{}

ant::TFSum::TFSum(const std::string& name, std::list<TF1*> fs, double xmin, double xmax):
    functions(fs),
    sum(BuildTF1(name, xmin,xmax))
{
    SetRange(xmin, xmax);
    SyncToSum();
}

ant::TFSum::~TFSum()
{

}

double ant::TFSum::operator()(const double *x, const double *p) const
{
    double s = 0.0;
    unsigned p_offset = 0;
    for(const auto& f : functions) {
        s += f->EvalPar(x, std::addressof(p[p_offset]));
        p_offset += f->GetNpar();
    }
    return s;
}

unsigned ant::TFSum::GetNpar() const { return sum->GetNpar(); }

void ant::TFSum::SetRange(double xmin, double xmax)
{
    for(const auto& f : functions) {
        f->SetRange(xmin,xmax);
    }

    sum->SetRange(xmin, xmax);
}

void ant::TFSum::GetRange(double &xmin, double &xmax) const
{
    sum->GetRange(xmin, xmax);
}

void ant::TFSum::SetNpx(int n)
{
    for(const auto& f : functions) {
        f->SetNpx(n);
    }

    sum->SetNpx(n);
}
