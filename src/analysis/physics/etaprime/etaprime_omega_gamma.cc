#include "etaprime_omega_gamma.h"
#include "plot/root_draw.h"
#include "utils/particle_tools.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include <cassert>
#include <numeric>

using namespace std;
using namespace ant::analysis;
using namespace ant::analysis::physics;




EtapOmegaG::EtapOmegaG(PhysOptPtr opts) : Physics("EtapOmegaG", opts)
{
    sig_steps = HistFac.makeTH1D("steps", "", "%", BinSettings(10),"steps");
}

EtapOmegaG::sig_perDecayHists_t::sig_perDecayHists_t(SmartHistFactory& HistFac_parent, const string& decaystring)
{
    auto directory_name = utils::ParticleTools::SanitizeDecayString(decaystring);
    SmartHistFactory HistFac(directory_name, HistFac_parent);
    HistFac.SetTitlePrefix(decaystring);

    BinSettings bins_im(1200);

    gggg = HistFac.makeTH1D("4#gamma IM","4#gamma IM / MeV","events",bins_im,"gggg");
    ggg = HistFac.makeTH1D("3#gamma IM","3#gamma IM / MeV","events",bins_im,"ggg");
    gg = HistFac.makeTH1D("2#gamma IM","2#gamma IM / MeV","events",bins_im,"gg");

    IM_etap_omega = HistFac.makeTH2D("#eta' vs. #omega IM",
                               "#eta' IM / MeV",
                               "#omega IM / MeV",
                               BinSettings(400, 600, 1200),
                               BinSettings(400, 400, 1000),
                               "IM_etap_omega");
    IM_pi0 = HistFac.makeTH1D("#pi^{0}","#pi^{0} IM / MeV","events",BinSettings(300, 0, 400),"IM_pi0");

    BinSettings bins_mm(300, 600, 1300);
    MM_gggg = HistFac.makeTH1D("M_{miss}","M_{miss} / MeV","events",bins_mm,"MM_gggg");
    MM_etap = HistFac.makeTH1D("M_{miss}^{cut}","M_{miss} / MeV","events",bins_mm,"MM_etap");

    Chi2_All = HistFac.makeTH1D("#chi^{2} all","#chi^{2}","",BinSettings(300,0,100),"Chi2_All");
    Chi2_Best = HistFac.makeTH1D("#chi^{2} minimal","#chi^{2}","",BinSettings(300,0,100),"Chi2_Min");

    Proton_ThetaPhi = HistFac.makeTH2D("p #delta(#theta-#phi)",
                                       "#delta#theta / degree",
                                       "#delta#phi / degree",
                                       BinSettings(100, -10, 10),
                                       BinSettings(100, -30, 30),
                                       "Proton_ThetaPhi"
                                       );
    Proton_Energy = HistFac.makeTH1D("p #delta(E)","#deltaE / MeV","",BinSettings(400,-50,350),"Proton_Energy");
}



void EtapOmegaG::ProcessEvent(const data::Event& event)
{
    const auto& data = event.Reconstructed();


    ProcessSig(getHistogram(event.MCTrue().ParticleTree(), sig_perDecayHists), data);
}

void EtapOmegaG::ProcessSig(sig_perDecayHists_t h,
                            const data::Event::Data& data)
{
    const auto nParticles = data.Particles().GetAll().size();

    sig_steps->Fill("Seen",1);

    if(nParticles != 5)
        return;

    sig_steps->Fill("nParticles==5",1);

    const auto& photons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles().Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 4)
        return;
    sig_steps->Fill("nPhotons==4",1);

    if(nProtons != 1)
        return;
    sig_steps->Fill("nProtons==1",1);
    const data::ParticlePtr& proton = protons.front();

    if(data.TriggerInfos().CBEenergySum() < 550)
        return;
    sig_steps->Fill("CBESum>550MeV",1);

    // gamma combinatorics
    assert(photons.size() == 4);
    utils::ParticleTools::FillIMCombinations(h.gg,   2, photons);
    utils::ParticleTools::FillIMCombinations(h.ggg,  3, photons);

    TLorentzVector photon_sum(0,0,0,0);
    for(const auto& p : photons) {
        photon_sum += *p;
    }
    h.gggg->Fill(photon_sum.M());

    // use tagged photon
    for(const auto& th : data.TaggerHits()) {
        /// \todo make prompt/random cut
        const TLorentzVector beam_target = th->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector mm = beam_target - photon_sum;
        h.MM_gggg->Fill(mm.M());
    }

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

    // photon assignment successful?
    if(result.Chi2>10)
        return;
    sig_steps->Fill("MinChi2<10",1);

    h.IM_etap_omega->Fill(result.EtaPrime.M(), result.Omega.M());
    h.IM_pi0->Fill(result.Pi0.M());
    h.Chi2_Best->Fill(result.Chi2);

    // use tagged photon
    for(const auto& th : data.TaggerHits()) {
        /// \todo make prompt/random cut
        const TLorentzVector beam_target = th->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector mm = beam_target - result.EtaPrime;
        h.MM_etap->Fill(mm.M());

        // don't assume we always have a proton...
        if(proton) {
            const double diff_Theta = mm.Theta() - proton->Theta();
            const double diff_Phi = mm.Phi() - proton->Phi();
            h.Proton_ThetaPhi->Fill(std_ext::radian_to_degree(diff_Theta), std_ext::radian_to_degree(diff_Phi));
            h.Proton_Energy->Fill(mm.E() - proton->E());
        }
    }
}


void EtapOmegaG::Finish()
{
    //steps->Scale(100.0/steps->GetBinContent(1));
}

void EtapOmegaG::ShowResult()
{
    for(const auto& it_map : sig_perDecayHists) {
        const sig_perDecayHists_t& h = it_map.second;
        if(h.IM_etap_omega->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it_map.first);
        c << sig_steps
          << h.gg << h.ggg << h.gggg << h.MM_gggg
          << h.Chi2_All << h.Chi2_Best
          << h.IM_pi0 << drawoption("colz") << h.IM_etap_omega
          << h.MM_etap
          << drawoption("colz") <<  h.Proton_ThetaPhi << h.Proton_Energy
          << endc;
    }
}






AUTO_REGISTER_PHYSICS(EtapOmegaG, "EtapOmegaG")
