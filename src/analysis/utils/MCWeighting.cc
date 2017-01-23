#include "MCWeighting.h"

#include "base/std_ext/string.h"

// use ROOT's function GSL wrapper...
#include "Math/SpecFuncMathMore.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

MCWeighting::MCWeighting() :
    database(MakeDatabase())
{

}

double MCWeighting::GetWeight(const TParticleTree_t& tree) const
{
    if(tree->Get()->Type() != ParticleTypeDatabase::BeamTarget)
        throw Exception("Root of ParticleTree must be beam particle");
    if(tree->Daughters().size() != 2)
        throw Exception("ParticleTree must be single particle production");

    // find the meson, so not the nucleon
    auto& meson = tree->Daughters().front()->Get()->Type() != ParticleTypeDatabase::Nucleon ?
                  tree->Daughters().front()->Get() : tree->Daughters().back()->Get();

    const auto it_d = database.find(addressof(meson->Type()));
    if(it_d == database.end())
        throw Exception(std_ext::formatter() <<
                        "Cannot find database entry for meson type " << meson->Type());

    auto& entry = it_d->second;

    // composite beamtarget particle, which is "incoming photon"+{proton,neutron,...}
    const double beamE = tree->Get()->Ek();


    // database coefficients are sorted in ascending BeamE
    // find the first item which has larger beamE
    using it_t = vector<coefficients_t>::const_iterator;
    it_t it_coeff_hi = std::find_if(entry.begin(), entry.end(),
                                 [beamE] (const coefficients_t& c) { return c.BeamE > beamE; });

    // handle corner cases of extrapolation
    if(it_coeff_hi == entry.begin()) {
        // beamE smaller than all entries
        it_coeff_hi = std::next(entry.begin());
    }
    else if(it_coeff_hi == entry.end()) {
        // beamE larger than all entries
        it_coeff_hi = std::prev(it_coeff_hi);
    }

    // get cosTheta in CM frame (rest frame of composite particle)
    LorentzVec meson_CM = *meson;
    meson_CM.Boost(-tree->Get()->BoostVector());
    const double cosTheta = std::cos(meson_CM.Theta());

    // get linearly interpolated legendre coefficients
    auto N = std_ext::NaN;
    {
        auto it_coeff_lo = std::prev(it_coeff_hi);

        // calculate N_lo / N_hi
        auto get_N = [cosTheta] (const it_t& it_coeff) {
            auto& coeffs = it_coeff->LegendreCoefficients;
            double sum = 0;
            for(auto l=0u;l<coeffs.size();l++)
                sum += coeffs[l]*ROOT::Math::legendre(l, cosTheta);
            return sum;
        };

        const double beamE_lo = it_coeff_lo->BeamE;
        const double beamE_hi = it_coeff_hi->BeamE;
        const double N_lo = get_N(it_coeff_lo);
        const double N_hi = get_N(it_coeff_hi);

        // slope
        const double m = (N_hi - N_lo)/(beamE_hi - beamE_lo);

        // calc interpolated N
        N = N_lo + m*(beamE - beamE_lo);
    }

    /// \todo this is not a weight, but finding the normalization is non-trivial
    return N;
}

MCWeighting::database_t MCWeighting::MakeDatabase()
{
    database_t d;

    // data for the EtaPrime was copied from P.Adlarson code
    // https://github.com/padlarson/a2GoAT/blob/AdlarsonAnalysis/configfiles/data.MC/etaprime_Legendrecoeff_effcorr_eta2g.txt
    // calculated by V.Kashevarov based on the
    // paper by S.Prakhov https://arxiv.org/abs/1701.04809
    d[addressof(ParticleTypeDatabase::EtaPrime)] = {
        {1450.2, {300.382,   29.304,  -14.001,   -4.470, -17.092}},
        {1456.8, {434.115,   68.364,   44.920,   -9.691, -19.235}},
        {1463.2, {496.643,   54.420,   52.968,  -12.282,  -4.047}},
        {1469.8, {549.032,   23.331,   53.547,   11.748,   9.741}},
        {1479.5, {1277.640, 111.258,   49.390,   -9.034, -49.627}},
        {1492.5, {1329.360, 131.958,  -23.574,    8.443, -92.103}},
        {1505.5, {1466.640, 203.069,  -57.099,  -31.691, -69.441}},
        {1518.5, {1531.590, 204.971, -132.569,  -19.244,   6.864}},
        {1531.5, {1587.570, 220.300,  -94.350,  -68.671,  17.392}},
        {1544.5, {1596.180, 261.441, -202.653,  -72.122, -62.119}},
        {1557.5, {1568.940, 271.070, -229.395,  -94.934,  13.390}},
        {1570.5, {1434.450, 244.342, -291.041, -130.754,  10.56}},
    };

    for(auto& i : d) {
        if(i.second.size()<2)
            throw Exception("Database coefficients must have at least two energy bins for linear interpolation");
        // ensure all coefficients have same size (useful for interpolation later)
        {
            auto it = i.second.begin();
            auto firstsize = it->LegendreCoefficients.size();
            if(firstsize==0)
                throw Exception("Found database entry with empty legendre coefficients");
            ++it;
            while(it != i.second.end()) {
                if(firstsize != it->LegendreCoefficients.size())
                    throw Exception("Found database entry with mismatching number of legendre coefficients");
                ++it;
            }
        }
        std::sort(i.second.begin(), i.second.end());
    }

    return d;
}
