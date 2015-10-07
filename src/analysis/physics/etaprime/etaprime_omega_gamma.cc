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

    gggg = HistFac.makeTH1D(decaystring+": 4#gamma IM","4#gamma IM / MeV","events",bins_im,pref+"_gggg");
    ggg = HistFac.makeTH1D(decaystring+": 3#gamma IM","3#gamma IM / MeV","events",bins_im,pref+"_ggg");
    gg = HistFac.makeTH1D(decaystring+": 2#gamma IM","2#gamma IM / MeV","events",bins_im,pref+"_gg");

    IM_etap_omega = HistFac.makeTH2D(decaystring+": #eta' vs. #omega IM",
                               "#eta' IM / MeV",
                               "#omega IM / MeV",
                               BinSettings(400, 600, 1200),
                               BinSettings(400, 400, 1000),
                               pref+"_IM_etap_omega");
    IM_pi0 = HistFac.makeTH1D(decaystring+": #pi^{0}","#pi^{0} IM / MeV","events",BinSettings(300, 0, 400),pref+"_IM_pi0");

    Chi2_All = HistFac.makeTH1D(decaystring+": #chi^{2} all","#chi^{2}","",BinSettings(300,0,100),pref+"_Chi2_All");
    Chi2_Best = HistFac.makeTH1D(decaystring+": #chi^{2} minimal","#chi^{2}","",BinSettings(300,0,100),pref+"_Chi2_Min");
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
        const auto Chi2_Pi0 =  std_ext::sqr((tmp.Pi0.M() - 126) / 15);
        const auto Chi2_Omega = std_ext::sqr((tmp.Omega.M() - 735) / 32);
        const auto Chi2_EtaPrime = std_ext::sqr((tmp.EtaPrime.M() - 895) / 27);

        tmp.Chi2 = Chi2_Pi0 + Chi2_Omega + Chi2_EtaPrime;

        h.Chi2_All->Fill(tmp.Chi2);

        if(tmp.Chi2<result.Chi2)
            result = tmp;
    }
    while(next_permutation(indices.begin(), indices.end(), comparer));

    if(result.Chi2<10) {
        h.IM_etap_omega->Fill(result.EtaPrime.M(), result.Omega.M());
        h.IM_pi0->Fill(result.Pi0.M());
        h.Chi2_Best->Fill(result.Chi2);
    }
}

void EtapOmegaG::ShowResult()
{
    for(const auto& it_map : perDecayHists) {
        const perDecayHists_t& h = it_map.second;
        if(h.IM_etap_omega->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it_map.first);
        c << h.gg << h.ggg << h.gggg
          << h.Chi2_All << h.Chi2_Best
          << h.IM_pi0 << drawoption("colz") << h.IM_etap_omega
          << endc;
    }
}






AUTO_REGISTER_PHYSICS(EtapOmegaG, "EtapOmegaG")
