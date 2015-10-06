#include "etaprime_omega_gamma.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include <cassert>

using namespace std;
using namespace ant::analysis;
using namespace ant::analysis::physics;


EtapOmegaG::EtapOmegaG(PhysOptPtr opts) : Physics("EtapOmegaG", opts)
{
    gggg = HistFac.makeTH1D("gggg","4#gamma IM / MeV","events",bins_im,"gggg");
    ggg = HistFac.makeTH1D("ggg","3#gamma IM / MeV","events",bins_im,"ggg");
    gg = HistFac.makeTH1D("gg","2#gamma IM / MeV","events",bins_im,"gg");

    IM_etap = HistFac.makeTH1D("EtaPrime","EtaPrime IM / MeV","events",bins_im,"IM_etap");
    IM_omega = HistFac.makeTH1D("Omega","Omega IM / MeV","events",bins_im,"IM_omega");
    IM_pi0 = HistFac.makeTH1D("Pi0","Pi0 IM / MeV","events",bins_im,"IM_pi0");

    Chi2_All = HistFac.makeTH1D("Chi2 all combinations","Chi2","",BinSettings(300,0,100),"Chi2_All");
    Chi2_Best = HistFac.makeTH1D("Chi2 best","Chi2","",BinSettings(300,0,100),"Chi2_Min");

    Chi2_Single_All = HistFac.makeTH2D("Chi2 components all combinations","Name","Chi2",
                                   BinSettings(3),
                                   BinSettings(300,0,100),
                                   "Chi2_Single_All");
    Chi2_Single_Best = HistFac.makeTH2D("Chi2 components best","Name","Chi2",
                                   BinSettings(3),
                                   BinSettings(300,0,100),
                                   "Chi2_Single_Best");
}

void EtapOmegaG::ProcessEvent(const data::Event& event)
{
    const auto& data = event.Reconstructed();

    const auto nParticles = data.Particles().GetAll().size();
    if(nParticles != 5)
        return;

    const auto& photons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles().Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 4 || nProtons != 1)
        return;

    // gamma combinatorics
    auto fill_combinations = [] (TH1* h, unsigned multiplicity, const data::ParticleList& particles) {
        for( auto comb = utils::makeCombination(particles,multiplicity); !comb.Done(); ++comb) {
             TLorentzVector sum(0,0,0,0);
             for(const auto& p : comb) {
                 sum += *p;
             }
             h->Fill(sum.M());
        }
    };

    fill_combinations(gg,   2, photons);
    fill_combinations(ggg,  3, photons);
    fill_combinations(gggg, 4, photons);

    // bottom-up assignment


    struct result_t {
        double Chi2 = std::numeric_limits<double>::infinity();
        double Chi2_Pi0;
        double Chi2_Omega;
        double Chi2_EtaPrime;
        data::ParticlePtr g_pi0_0;
        data::ParticlePtr g_pi0_1;
        data::ParticlePtr g_omega;
        data::ParticlePtr g_etap;
        TLorentzVector Pi0;
        TLorentzVector Omega;
        TLorentzVector EtaPrime;
    };

    result_t result; // best chi2

    assert(photons.size() == 4);
    vector<unsigned> indices = {0,1,2,3};
    auto comparer = [] (unsigned a, unsigned b) {
        // make 0 and 1 look equal
        if(a==0 && b==1)
            return false;
        return a<b;
    };

    do {
        VLOG(9) << indices;
        // the i vector tells us what particle
        // should be used as daughter particle
        // 0,1 : from pi0
        // 2   : from omega
        // 3   : from eta'

        result_t tmp;

        for(unsigned i=0;i<indices.size();i++) {
            switch(indices[i]) {
            case 0: tmp.g_pi0_0 = photons[i]; break;
            case 1: tmp.g_pi0_1 = photons[i]; break;
            case 2: tmp.g_omega = photons[i]; break;
            case 3: tmp.g_etap  = photons[i]; break;
            }
        }

        tmp.Pi0 = *(tmp.g_pi0_0) + *(tmp.g_pi0_1);

        tmp.Omega = tmp.Pi0 + *tmp.g_omega;
        tmp.EtaPrime = tmp.Omega + *tmp.g_etap;

        // means/sigma extracted from gg/ggg/gggg histograms
        tmp.Chi2_Pi0 =  std_ext::sqr((tmp.Pi0.M() - 126) / 15);
        tmp.Chi2_Omega = std_ext::sqr((tmp.Omega.M() - 735) / 32);
        tmp.Chi2_EtaPrime = std_ext::sqr((tmp.EtaPrime.M() - 895) / 27);

        Chi2_Single_All->Fill("Pi0",tmp.Chi2_Pi0, 1.0);
        Chi2_Single_All->Fill("Omega",tmp.Chi2_Omega, 1.0);
        Chi2_Single_All->Fill("EtaP",tmp.Chi2_EtaPrime, 1.0);


        tmp.Chi2 = tmp.Chi2_Pi0 + tmp.Chi2_Omega + tmp.Chi2_EtaPrime;

        Chi2_All->Fill(tmp.Chi2);

        if(tmp.Chi2<result.Chi2)
            result = tmp;
    }
    while(next_permutation(indices.begin(), indices.end(), comparer));

    IM_etap->Fill(result.EtaPrime.M());
    IM_omega->Fill(result.Omega.M());
    IM_pi0->Fill(result.Pi0.M());
    Chi2_Best->Fill(result.Chi2);
    Chi2_Single_Best->Fill("Pi0",result.Chi2_Pi0, 1.0);
    Chi2_Single_Best->Fill("Omega",result.Chi2_Omega, 1.0);
    Chi2_Single_Best->Fill("EtaP",result.Chi2_EtaPrime, 1.0);



}

void EtapOmegaG::ShowResult()
{
    canvas(GetName()) << gg << ggg << gggg
                      << IM_pi0 << IM_omega << IM_etap
                      << Chi2_All << Chi2_Best
                      << drawoption("colz") << Chi2_Single_All << Chi2_Single_Best
                      << endc;
}


AUTO_REGISTER_PHYSICS(EtapOmegaG, "EtapOmegaG")
