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

// from http://dx.doi.org/10.1140/epja/i2003-10061-y
// does not include tagger spectrum
const MCWeighting::item_t MCWeighting::Omega = {
    ParticleTypeDatabase::Omega,
    MCWeighting::SanitizeDatabase({
        {{1144.2553238418518,	1174.999205771379}, 	{0.491519228648931, 0.220528891613139, 0.023578457249135, 0.0144933665460678, 0.00448044514746282, -0.0603381155962589}},
        {{1174.999205771379,    1206.0332487487183},	{0.582823716700338, 0.282032109536378, -0.00734705360199404, 0.132650715871278, 0.0283473481205772, 0.0506134608146485}},
        {{1206.0332487487183,	1249.7769584488785},	{0.633658642529541, 0.288433523620573, 0.097987927550612, 0.151117989066987, 0.117464702641173, 0.031715570901563}},
        {{1249.7769584488785,	1299.9052769898751},	{0.692304441611994, 0.380053144915676, 0.200434782117419, 0.183199466994578, 0.13250836462011, -0.00712591609835936}},
        {{1299.9052769898751,	1349.7694130892687},	{0.624483973761887, 0.401298830837046, 0.210593477609198, 0.277966179123474, 0.169902308090458, 0.0282337638018723}},
        {{1349.7694130892687,	1399.3286003188534},	{0.580254967503796, 0.429851112483861, 0.306674611328718, 0.336544298107312, 0.255332277421272, 0.191441376939843}},
        {{1399.3286003188534,	1437.439082403918},	    {0.545898136607653, 0.458586181164415, 0.412133267818453, 0.531182344272846, 0.271929650742513, 0.143559918866463}},
        {{1437.439082403918,	1462.7215951901937},	{0.486308861698685, 0.480966752771833, 0.405723864173806, 0.291563818784434, 0.154323042405376, 0.00308279619835812}},
        {{1462.7215951901937,	1499.4209726754311},	{0.497135376708792, 0.461983610603222, 0.443731579123435, 0.572872101105675, 0.350003488453786, 0.181533842729711}},
        {{1499.4209726754311,	1549.9277795253497},	{0.514079539437838, 0.584611713919711, 0.578462551783947, 0.674652063739541, 0.375789808688014, 0.148662793197016}},
        {{1549.9277795253497,	1600.0241243971652},	{0.483419245891268, 0.616851705001912, 0.65752616160673, 0.707671280941172, 0.511081861164857, 0.216591998887943}},
        {{1600.0241243971652,	1649.6716388878608},	{0.58000830141939, 0.839313691646204, 1.05040985804401, 0.928172448380859, 0.721133476371126, 0.23839347563335}},
        {{1649.6716388878608,	1699.9077353387183},	{0.548659891440818, 0.852265435849442, 0.962189377923219, 0.893180741559411, 0.584274978669801, 0.325875737671703}},

    })
};

// only valid in EPT range
// from https://arxiv.org/abs/1506.08849
// does not include tagger spectrum
const MCWeighting::item_t MCWeighting::Pi0 = {
    ParticleTypeDatabase::Pi0,
    {
        {{1400,1500}, 	{1.04444,0.471311,1.39205,0.827052,1.34353,0.203569}},
        {{1500,1600}, 	{1.04444,0.471311,1.39205,0.827052,1.34353,0.203569}},
    }
};

// from https://arxiv.org/abs/1701.04809
// does not include tagger spectrum
const MCWeighting::item_t MCWeighting::Eta = {
    ParticleTypeDatabase::Omega,
    MCWeighting::SanitizeDatabase({
        {{1402.32,	1452.59},	{0.150486582529627,	0.0891658300725101	,-0.0260755674728463 ,   -0.0198294310291532,	0.0104289851916939,	 -0.0174449879242664,	-0.024565612765031,  -0.0160543372119484}},
        {{1452.59,	1503.52},	{0.150015817651979,	0.0891619718173603	,-0.0248942512155738 ,   -0.016829583701192,	0.0120791150218525,  -0.0153853326533159,	-0.0231266413402218, -0.0154326183015057}},
        {{1503.52,	1550.97},	{0.11400430990043 , 0.10367371804368	,0.00950627895434919 ,   -0.00570912336165498, -0.00688196192920969, -0.0294076971497731,	-0.0421575979497328, -0.0285574279526046}},
        {{1550.97,	1588.49},	{0.114237939359211,	0.112707134960017	,0.0266870613695321	 ,    0.0144464343986868,	0.00465993319736669, -0.0175052336164756,	-0.0376896812319864, -0.0211751434157216}},
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








