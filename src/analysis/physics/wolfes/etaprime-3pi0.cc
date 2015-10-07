#include "etaprime-3pi0.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "base/std_ext/math.h"

#include <algorithm>


using namespace std;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Etap3pi0::Etap3pi0(PhysOptPtr opts) : Physics("EtapOmegaG", opts)
{
    BinSettings bs = BinSettings(1600);
    hNgammaMC  = HistFac.makeTH1D("# gamma MC true","# #gamma","# events",BinSettings(8));
    hNgamma    = HistFac.makeTH1D("# gamma","# #gamma","# events",BinSettings(8));
    h2g        = HistFac.makeTH1D("2 #gamma","2#gamma IM [MeV]","#",bs,"gg");
    h6g        = HistFac.makeTH1D("6 #gamma","6#gamma IM [MeV]","#",bs,"gggggg");

    IM_etap    = HistFac.makeTH1D("EtaPrime","EtaPrime IM [MeV]","events",bs,"IM_etap");
    IM_pi0     = HistFac.makeTH1D("Pi0","Pi0 IM [MeV]","events",bs,"IM_pi0");

    IM_vs_chi2 = HistFac.makeTH2D("#chi2","Pi0 candidate mass [MeV]","#chi2",bs,BinSettings(100));

}

void Etap3pi0::ProcessEvent(const data::Event& event)
{
    const auto& data   = event.Reconstructed();
    const auto& mcdata = event.MCTrue();


    const auto& photons   = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& mcphotons = mcdata.Particles().Get(ParticleTypeDatabase::Photon);

    hNgamma->Fill(photons.size());
    hNgammaMC->Fill(mcphotons.size());

    if ( photons.size() != 6)
        return;

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

    struct pi0candidate
    {
        double chi2;
        TLorentzVector V;
    };
    vector<pi0candidate> pi0List;


    for ( auto g1 = 0 ; g1 < 5 ; ++g1)
        for ( auto g2 = g1 + 1; g2 < 6 ; ++g2)
        {
            pi0candidate pc;
            pc.V = *(photons.at(g1)) + *(photons.at(g2));
            pc.chi2 = std_ext::sqr(pc.V.M() - ParticleTypeDatabase::Pi0.Mass() );
            pi0List.push_back(pc);
        }

    std::sort(pi0List.begin(),pi0List.end(),
              [](const pi0candidate& a,
                 const pi0candidate& b) -> bool
                 {
                     return a.chi2 < b.chi2;
                 } );

    for (const auto& pcan: pi0List)
        IM_vs_chi2->Fill(pcan.V.M(),pcan.chi2);



    if (pi0List.size() < 3) return;
    TLorentzVector etapV(0,0,0,0);
    for (auto i = 0 ; i < 3 ; ++i)
    {
        IM_pi0->Fill(pi0List.at(i).V.M());
        etapV += pi0List.at(i).V;
    }
    IM_etap->Fill(etapV.M());

}


void Etap3pi0::ShowResult()
{
    canvas(GetName()) << h2g << h6g
                      << IM_pi0 << IM_etap
                      << hNgamma << hNgammaMC
                      << endc;
}


AUTO_REGISTER_PHYSICS(Etap3pi0, "Etap3pi0")
