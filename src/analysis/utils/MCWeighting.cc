#include "MCWeighting.h"

#include "utils/ParticleTools.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"


// use ROOT's function GSL wrapper...
#include "Math/SpecFuncMathMore.h"

#include <numeric>

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

const string MCWeighting::treeName = "MCWeighting";

// data for the EtaPrime was copied from P.Adlarson code
// but actually merged from provided files by mail (uses Sergey's original binning not Viktors bin center positions
// https://github.com/padlarson/a2GoAT/blob/AdlarsonAnalysis/configfiles/data.MC/etaprime_Legendrecoeff_effcorr_eta2g.txt
// calculated by V.Kashevarov based on the
// paper by S.Prakhov https://arxiv.org/abs/1701.04809
const MCWeighting::item_t MCWeighting::EtaPrime = {
    ParticleTypeDatabase::EtaPrime,
    MCWeighting::SanitizeDatabase({
        {{1447.0, 1453.5}, {300.382,   29.304,  -14.001,   -4.470, -17.092}},
        {{1453.5, 1460.0}, {434.115,   68.364,   44.920,   -9.691, -19.235}},
        {{1460.0, 1466.5}, {496.643,   54.420,   52.968,  -12.282,  -4.047}},
        {{1466.5, 1473.0}, {549.032,   23.331,   53.547,   11.748,   9.741}},
        {{1473.0, 1486.0}, {1277.640, 111.258,   49.390,   -9.034, -49.627}},
        {{1486.0, 1499.0}, {1329.360, 131.958,  -23.574,    8.443, -92.103}},
        {{1499.0, 1512.0}, {1466.640, 203.069,  -57.099,  -31.691, -69.441}},
        {{1512.0, 1525.0}, {1531.590, 204.971, -132.569,  -19.244,   6.864}},
        {{1525.0, 1538.0}, {1587.570, 220.300,  -94.350,  -68.671,  17.392}},
        {{1538.0, 1551.0}, {1596.180, 261.441, -202.653,  -72.122, -62.119}},
        {{1551.0, 1564.0}, {1568.940, 271.070, -229.395,  -94.934,  13.390}},
        {{1564.0, 1577.0}, {1434.450, 244.342, -291.041, -130.754,  10.56}},
    })
};

MCWeighting::MCWeighting(const HistogramFactory& histFac, const item_t& item) :
    Item(item),
    HistFac(histFac)
{
}

MCWeighting::database_t MCWeighting::SanitizeDatabase(database_t d)
{
    if(d.size()<2)
        throw Exception("Database coefficients must have at least two energy bins for linear interpolation");

    // just to be sure
    sort(d.begin(), d.end());

    // normalize Legendre coefficients to width of beamE bin
    // (as database might have differently sized energy bins...)
    for(auto& e : d)
        for(auto& c : e.LegendreCoefficients)
            c /= e.BeamE.Length();

    return d;
}

double MCWeighting::GetBeamE(const TParticleTree_t& tree)
{
    // composite beamtarget particle, which is "incoming photon"+{proton,neutron,...}
    return tree->Get()->Ek();
}

double MCWeighting::GetCosTheta(const TParticleTree_t& tree)
{
    // find the meson, so not the nucleon
    auto& meson = tree->Daughters().front()->Get()->Type() != ParticleTypeDatabase::Nucleon ?
                                                                  tree->Daughters().front()->Get() : tree->Daughters().back()->Get();
    // get cosTheta in CM frame (rest frame of composite particle)
    LorentzVec meson_CM = *meson;
    meson_CM.Boost(-tree->Get()->BoostVector());
    return std::cos(meson_CM.Theta());
}

double MCWeighting::GetN(const double beamE, const double cosTheta) const
{
    // database coefficients are sorted in ascending BeamE
    // find the first item which has larger beamE
    const auto& db = Item.Database;
    using it_t = vector<coefficients_t>::const_iterator;
    it_t it_coeff_hi = std::find_if(db.begin(), db.end(),
                                    [beamE] (const coefficients_t& c) { return c.BeamE.Center() > beamE; });

    // handle corner cases of extrapolation
    if(it_coeff_hi == db.begin()) {
        // beamE smaller than all entries
        it_coeff_hi = std::next(db.begin());
    }
    else if(it_coeff_hi == db.end()) {
        // beamE larger than all entries
        it_coeff_hi = std::prev(it_coeff_hi);
    }

    // get linearly interpolated legendre coefficients
    auto it_coeff_lo = std::prev(it_coeff_hi);

    // calculate N_lo / N_hi
    auto get_N = [cosTheta] (const it_t& it_coeff) {
        auto& coeffs = it_coeff->LegendreCoefficients;
        double sum = 0;
        for(auto l=0u;l<coeffs.size();l++)
            sum += coeffs[l]*ROOT::Math::legendre(l, cosTheta);
        return sum;
    };

    const double beamE_lo = it_coeff_lo->BeamE.Center();
    const double beamE_hi = it_coeff_hi->BeamE.Center();
    const double N_lo = get_N(it_coeff_lo);
    const double N_hi = get_N(it_coeff_hi);

    // slope
    const double m = (N_hi - N_lo)/(beamE_hi - beamE_lo);

    // calc energy-interpolated N
    return N_lo + m*(beamE - beamE_lo);
}

void MCWeighting::SetParticleTree(const TParticleTree_t& tree)
{
    if(tree && tree->Get()->Type() != ParticleTypeDatabase::BeamTarget)
        throw Exception("Root of ParticleTree must be beam particle");

    // lazy init of tree
    if(t.Tree == nullptr)
        t.CreateBranches(HistFac.makeTTree(treeName+"_UNFINISHED"));

    // check if it the specified meson was produced
    if(ParticleTools::FindParticle(Item.Type, tree, 1) &&
       tree->Daughters().size() == 2) {
        last_N = GetN(GetBeamE(tree), GetCosTheta(tree));
        N_sum += last_N;
        nParticleTrees++;
    }
    else {
        last_N = std_ext::NaN;
    }
}

void MCWeighting::Fill()
{
    if(!t.Tree)
        return;

    t.MCWeight = last_N;
    t.Tree->Fill();
}

void MCWeighting::Finish()
{
    if(!t.Tree)
        return;

    tree_t t_norm;
    t_norm.CreateBranches(HistFac.makeTTree(treeName));

    for(auto entry = 0;entry<t.Tree->GetEntries(); entry++) {
        t.Tree->GetEntry(entry);
        t_norm.MCWeight = isfinite(t.MCWeight) ? (t.MCWeight * nParticleTrees)/N_sum : 1.0;
        t_norm.Tree->Fill();
    }

    // dispose temp tree
    delete t.Tree;
    t.Tree = nullptr;

    // remember tree for possible friending
    treeWeighted = t_norm.Tree;
}

bool MCWeighting::FriendTTree(TTree* tree)
{
    if(treeWeighted) {
        tree->AddFriend(treeWeighted, "", kTRUE);
        return true;
    }
    return false;
}








