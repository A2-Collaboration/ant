#include "Math.h"
#include "base/math_functions/AsymGaus.h"
#include "base/math_functions/CrystalBall.h"

#include "TH1.h"

using namespace ant;
using namespace std;

TF1* Math::AsymGaus()
{
    return math::AsymGaus::GetTF1();
}

TF1* Math::CrystalBall()
{
    return math::CrystalBall::GetTF1();
}

void TFSum::syncTF1Par(TF1* src, TF1* dst) {
    if(src->GetNpar() != dst->GetNpar())
        throw runtime_error("Not equal parameters");
    for(int j=0;j<src->GetNpar();j++)
        syncTF1Par(src, dst, j, j);
}

void TFSum::syncTF1Par(TF1* src, TF1* dst, int src_j, int dst_j) {
    dst->SetParameter(dst_j, src->GetParameter(src_j));
    double low, high;
    src->GetParLimits(src_j, low, high);
    dst->SetParLimits(dst_j, low, high);
    dst->SetParError(dst_j, src->GetParError(src_j));
    dst->SetParName(dst_j,  src->GetParName(src_j));
}

void TFSum::SyncToSum()
{
    int par = 0;
    for(const auto& f : functions) {
        for(int j=0; j<f->GetNpar(); ++j) {
            syncTF1Par(f, sum, j, par);
            ++par;
        }
    }
}

void TFSum::SyncToFcts()
{
    int par = 0;
    for(auto& f : functions) {
        for(int j=0; j<f->GetNpar(); ++j) {
            syncTF1Par(sum, f, par, j);
            ++par;
        }
    }
}

void TFSum::Draw()
{
    sum->Draw("same");
    for(auto& f : functions) {
        f->Draw("same");
    }
}

void TFSum::Draw(const string& option) const {
    (void)option;
    sum->Draw("same");
    for(auto& f : functions) {
        f->Draw("same");
    }
}

TF1* TFSum::BuildTF1(const string& name, const double min, const double max) const
{
    unsigned pars = 0;
    for(const auto& f : functions) {
        pars += f->GetNpar();
    }

    return new TF1(name.c_str(), this, min, max, pars);
}

TFSum::TFSum(const string &name, TF1 *f1, double xmin, double xmax):
    TFSum(name, list<TF1*>({f1}), xmin, xmax)
{}

TFSum::TFSum(const string &name, TF1 *f1, TF1 *f2, double xmin, double xmax):
    TFSum(name, {f1,f2}, xmin, xmax)
{}

TFSum::TFSum(const string &name, TF1 *f1, TF1 *f2, TF1* f3, double xmin, double xmax):
    TFSum(name, {f1,f2,f3}, xmin, xmax)
{}

TFSum::TFSum(const string& name, list<TF1*> fs, double xmin, double xmax):
    functions(fs),
    sum(BuildTF1(name, xmin,xmax))
{
    SetRange(xmin, xmax);
    SyncToSum();
}

TFSum::~TFSum()
{

}

double TFSum::operator()(const double *x, const double *p) const
{
    double s = 0.0;
    unsigned p_offset = 0;
    for(const auto& f : functions) {
        s += f->EvalPar(x, addressof(p[p_offset]));
        p_offset += f->GetNpar();
    }
    return s;
}

unsigned TFSum::GetNpar() const { return sum->GetNpar(); }

void TFSum::SetRange(double xmin, double xmax)
{
    for(const auto& f : functions) {
        f->SetRange(xmin,xmax);
    }

    sum->SetRange(xmin, xmax);
}

void TFSum::GetRange(double &xmin, double &xmax) const
{
    sum->GetRange(xmin, xmax);
}

void TFSum::SetNpx(int n)
{
    for(const auto& f : functions) {
        f->SetNpx(n);
    }

    sum->SetNpx(n);
}

TFitResultPtr TFSum::FitRanged(TH1* h, TF1* f, double x_low, double x_high, const string& fitopts)
{
    return FitRanged(h, f, {{x_low, x_high}}, fitopts);
}

TFitResultPtr TFSum::FitRanged(TH1* h, TF1* f,
                               double x1_low, double x1_high, double x2_low, double x2_high,
                               const string& fitopts)
{
    return FitRanged(h, f, {{x1_low, x1_high},{x2_low, x2_high}}, fitopts);
}

TF1* TFSum::MakeRanged(TF1* f, double x_low, double x_high)
{
    return MakeRanged(f, {{x_low, x_high}});
}

TF1* TFSum::MakeRanged(TF1* f, double x1_low, double x1_high, double x2_low, double x2_high)
{
    return MakeRanged(f, {{x1_low, x1_high},{x2_low, x2_high}});
}

TFitResultPtr TFSum::FitRanged(TH1* h, TF1* f, const PiecewiseInterval<double>& range, const string& fitopts)
{
    auto f_ranged = MakeRanged(f, range);
    auto r = h->Fit(f_ranged, fitopts.c_str());
    syncTF1Par(f_ranged, f);
    return r;
}

TF1* TFSum::MakeRanged(TF1* f, const PiecewiseInterval<double>& range)
{
    return MakeFiltered(f, [range] (double* x) {
        return range.Contains(x[0]);
    });
}
