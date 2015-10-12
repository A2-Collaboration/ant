#include "etaprime-3pi0.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "base/std_ext/math.h"
#include "data/Particle.h"

#include <algorithm>
#include <cassert>

#include "TTree.h"

using namespace std;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;



Etap3pi0::Etap3pi0(PhysOptPtr opts) :
    Physics("EtapOmegaG", opts),
      dataset(opts->GetOption("dataset"))
{
    BinSettings bs = BinSettings(1600);
    hNgammaMC  = HistFac.makeTH1D("# gamma MC true","# #gamma","# events",BinSettings(8));
    hNgamma    = HistFac.makeTH1D("# gamma","# #gamma","# events",BinSettings(8));
    h2g        = HistFac.makeTH1D("2 #gamma","2#gamma IM [MeV]","#",bs,"gg");
    h6g        = HistFac.makeTH1D("6 #gamma","6#gamma IM [MeV]","#",bs,"gggggg");

    IM_etap    = HistFac.makeTH1D("EtaPrime","EtaPrime IM [MeV]","events",bs,"IM_etap");
    IM_pi0     = HistFac.makeTH1D("Pi0","Pi0 IM [MeV]","events",bs,"IM_pi0");

    tree       = new TTree("data","data");

    for ( auto& pi0: p0cands)
        pi0.SetBranches(tree,(dataset + "pi0cands"));
    for ( auto& pi0: p0best)
        pi0.SetBranches(tree,(dataset + "pi0best"));

    MMproton.SetBranches(tree,(dataset + "MMproton"));
    etaprime.SetBranches(tree,(dataset + "etaprime"));

}

void Etap3pi0::ProcessEvent(const data::Event& event)
{
    const auto& data   = event.Reconstructed();
    const auto& mcdata = event.MCTrue();


    const auto& photons   = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& mcphotons = mcdata.Particles().Get(ParticleTypeDatabase::Photon);

    hNgamma->Fill(photons.size());
    hNgammaMC->Fill(mcphotons.size());

    auto fill_combinations = [] (TH1* h, unsigned multiplicity, const data::ParticleList& particles) {
        for( auto comb = utils::makeCombination(particles,multiplicity); !comb.Done(); ++comb) {
             TLorentzVector sum(0,0,0,0);
             for(const auto& p : comb) {
                 sum += *p;
             }
             h->Fill(sum.M());
        }
    };

    fill_combinations(h2g, 2, photons);
    fill_combinations(h6g, 6, photons);

    struct result_t {
        double Chi2 = std::numeric_limits<double>::infinity();
        vector<data::ParticlePtr> g_pi0;
        vector<TLorentzVector> Pi0;
        result_t() : g_pi0(6), Pi0(3) {}
    };

    assert(photons.size() == 6);
    vector<unsigned> indices = {0,1,2,3,4,5};
    auto comparer = [] (unsigned a, unsigned b) {
        // make 0/1 and 2/3 and 4/5 look equal
        if(a==0 && b==1)
            return false;
        if(a==2 && b==3)
            return false;
        if(a==4 && b==5)
            return false;
        return a<b;
    };

    result_t result; // best chi2


    do {
        // the indices vector tells us what particle
        // should be used as daughter particle
        // 0,1 : from first pi0
        // 2,3 : from second pi0
        // 4,5 : from third pi0
        result_t tmp;

        for(unsigned i=0;i<indices.size();i++) {
            tmp.g_pi0[indices[i]] = photons[i];
        }

        tmp.Chi2 = 0;
        for(unsigned i=0;i<tmp.Pi0.size();i++) {
            tmp.Pi0[i] = *(tmp.g_pi0[2*i]) + *(tmp.g_pi0[2*i+1]);
            tmp.Chi2 += std_ext::sqr((tmp.Pi0[i].M() - 126) / 15);
        }
        for (const auto& p: tmp.Pi0 )
        {
            IM_pi0->Fill(p.M());
        }

        if(tmp.Chi2<result.Chi2)
            result = move(tmp);
    }
    while(next_permutation(indices.begin(), indices.end(), comparer));

    TLorentzVector etap(0,0,0,0);
    for (const auto& p: result.Pi0 )
    {
        etap += p;
    }
    IM_etap->Fill(etap.M());




}


void Etap3pi0::ShowResult()
{
    canvas(GetName()) << h2g << h6g
                      << IM_pi0 << IM_etap
                      << hNgamma << hNgammaMC
                      << falseNgamma
                      << endc;
}

Etap3pi0::ParticleVars::ParticleVars(const TLorentzVector& lv, const ParticleTypeDatabase::Type& type) noexcept
{
    IM    = lv.M();
    Theta = radian_to_degree(lv.Theta());
    Phi   = radian_to_degree(lv.Phi());
    E     = lv.E() - type.Mass();
}

Etap3pi0::ParticleVars::ParticleVars(const Particle& p) noexcept
{
    IM    = p.M();
    Theta = radian_to_degree(p.Theta());
    Phi   = radian_to_degree(p.Phi());
    E     = p.Ek();
}

void Etap3pi0::ParticleVars::SetBranches(TTree* tree, const string& name)
{
    tree->Branch((name+"IM").c_str(), &IM);
    tree->Branch((name+"Theta").c_str(), &Theta);
    tree->Branch((name+"Phi").c_str(),&Phi);
    tree->Branch((name+"E").c_str(),  &E);
}

AUTO_REGISTER_PHYSICS(Etap3pi0, "Etap3pi0")
