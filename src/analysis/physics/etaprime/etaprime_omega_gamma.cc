#include "etaprime_omega_gamma.h"
#include "plot/root_draw.h"
#include "utils/particle_tools.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include <cassert>

using namespace std;
using namespace ant::analysis;
using namespace ant::analysis::physics;


EtapOmegaG::EtapOmegaG(PhysOptPtr opts) : Physics("EtapOmegaG", opts)
{

}

EtapOmegaG::perDecayHists_t::perDecayHists_t(SmartHistFactory& HistFac, const string& decaystring)
{
    auto pref = utils::ParticleTools::SanitizeDecayString(decaystring);

    BinSettings bins_im(1200);

    gggg = HistFac.makeTH1D(decaystring+": gggg","4#gamma IM / MeV","events",bins_im,pref+"_gggg");
    ggg = HistFac.makeTH1D(decaystring+": ggg","3#gamma IM / MeV","events",bins_im,pref+"_ggg");
    gg = HistFac.makeTH1D(decaystring+": gg","2#gamma IM / MeV","events",bins_im,pref+"_gg");

    IM_etap = HistFac.makeTH1D(decaystring+": EtaPrime","EtaPrime IM / MeV","events",bins_im,pref+"_IM_etap");
    IM_omega = HistFac.makeTH1D(decaystring+": Omega","Omega IM / MeV","events",bins_im,pref+"_IM_omega");
    IM_pi0 = HistFac.makeTH1D(decaystring+": Pi0","Pi0 IM / MeV","events",bins_im,pref+"_IM_pi0");

    Chi2_All = HistFac.makeTH1D(decaystring+": Chi2 all combinations","Chi2","",BinSettings(300,0,100),pref+"_Chi2_All");
    Chi2_Best = HistFac.makeTH1D(decaystring+": Chi2 best","Chi2","",BinSettings(300,0,100),pref+"_Chi2_Min");

    Chi2_Single_All = HistFac.makeTH2D(decaystring+": Chi2 components all combinations","Name","Chi2",
                                   BinSettings(3),
                                   BinSettings(300,0,100),
                                   pref+"_Chi2_Single_All");
    Chi2_Single_Best = HistFac.makeTH2D(decaystring+": Chi2 components best","Name","Chi2",
                                   BinSettings(3),
                                   BinSettings(300,0,100),
                                   pref+"_Chi2_Single_Best");
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



    const string& decaystring = utils::ParticleTools::GetDecayString(event.MCTrue().Intermediates().GetAll());

    // search map only once even on insert
    auto it_h = perDecayHists.lower_bound(decaystring);
    if(it_h == perDecayHists.end() || perDecayHists.key_comp()(decaystring, it_h->first)) {
        // decaystring does not exist
        it_h = perDecayHists.emplace_hint(it_h, decaystring, perDecayHists_t(HistFac, decaystring));
    }
    const perDecayHists_t& h = it_h->second;

    // gamma combinatorics
    utils::ParticleTools::FillIMCombinations(h.gg,   2, photons);
    utils::ParticleTools::FillIMCombinations(h.ggg,  3, photons);
    utils::ParticleTools::FillIMCombinations(h.gggg, 4, photons);

    // bottom-up assignment of photons

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

        h.Chi2_Single_All->Fill("Pi0",tmp.Chi2_Pi0, 1.0);
        h.Chi2_Single_All->Fill("Omega",tmp.Chi2_Omega, 1.0);
        h.Chi2_Single_All->Fill("EtaP",tmp.Chi2_EtaPrime, 1.0);


        tmp.Chi2 = tmp.Chi2_Pi0 + tmp.Chi2_Omega + tmp.Chi2_EtaPrime;

        h.Chi2_All->Fill(tmp.Chi2);

        if(tmp.Chi2<result.Chi2)
            result = tmp;
    }
    while(next_permutation(indices.begin(), indices.end(), comparer));

    if(result.Chi2<10) {
        h.IM_etap->Fill(result.EtaPrime.M());
        h.IM_omega->Fill(result.Omega.M());
        h.IM_pi0->Fill(result.Pi0.M());
        h.Chi2_Best->Fill(result.Chi2);
        h.Chi2_Single_Best->Fill("Pi0",result.Chi2_Pi0, 1.0);
        h.Chi2_Single_Best->Fill("Omega",result.Chi2_Omega, 1.0);
        h.Chi2_Single_Best->Fill("EtaP",result.Chi2_EtaPrime, 1.0);
    }
}

void EtapOmegaG::ShowResult()
{
    for(const auto& it_map : perDecayHists) {
        canvas c(GetName()+": "+it_map.first);
        const perDecayHists_t& h = it_map.second;
        c << h.gg << h.ggg << h.gggg
          << h.IM_pi0 << h.IM_omega << h.IM_etap
          << h.Chi2_All << h.Chi2_Best
          << drawoption("colz") << h.Chi2_Single_All << h.Chi2_Single_Best;
        c << endc;
    }
}






AUTO_REGISTER_PHYSICS(EtapOmegaG, "EtapOmegaG")
