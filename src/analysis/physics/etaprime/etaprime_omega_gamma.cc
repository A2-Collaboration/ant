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




EtapOmegaG::EtapOmegaG(PhysOptPtr opts) : Physics("EtapOmegaG", opts),
    HistFac_sig("Sig",HistFac),
    HistFac_ref("Ref",HistFac)
{
    sig_steps = HistFac_sig.makeTH1D("Steps: Signal channel", "", "#", BinSettings(10),"sig_steps");
    ref_steps = HistFac_ref.makeTH1D("Steps: Reference channel", "", "#", BinSettings(10),"ref_steps");

    treeSig = HistFac_sig.makeTTree("treeSig");
    treeRef = HistFac_ref.makeTTree("treeRef");


}

EtapOmegaG::sig_perDecayHists_t::sig_perDecayHists_t(SmartHistFactory& HistFac_parent,
                                                     const string& decaystring)
{
    auto directory_name = utils::ParticleTools::SanitizeDecayString(decaystring);
    SmartHistFactory HistFac(directory_name, HistFac_parent);
    HistFac.SetTitlePrefix(decaystring);

    BinSettings bins_im(1200);

    gggg = HistFac.makeTH1D("4#gamma IM","4#gamma IM / MeV","events",bins_im,"gggg");
    ggg = HistFac.makeTH1D("3#gamma IM","3#gamma IM / MeV","events",bins_im,"ggg");
    gg = HistFac.makeTH1D("2#gamma IM","2#gamma IM / MeV","events",bins_im,"gg");

    Proton_Copl = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");

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

EtapOmegaG::ref_perDecayHists_t::ref_perDecayHists_t(SmartHistFactory& HistFac_parent,
                                                     const string& decaystring)
{
    auto directory_name = utils::ParticleTools::SanitizeDecayString(decaystring);
    SmartHistFactory HistFac(directory_name, HistFac_parent);
    HistFac.SetTitlePrefix(decaystring);

    BinSettings bins_im(1200);

    gg = HistFac.makeTH1D("2#gamma IM","2#gamma IM / MeV","events",bins_im,"gg");

    Proton_Copl = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");

    IM_etap = HistFac.makeTH1D("#eta' IM","#eta' IM / MeV","events",BinSettings(300, 700, 1100),"IM_etap");

    BinSettings bins_mm(300, 600, 1300);
    MM_etap = HistFac.makeTH1D("M_{miss}","M_{miss} / MeV","events",bins_mm,"MM_etap");


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

    ProcessRef(event.MCTrue().ParticleTree(), data);
    ProcessSig(event.MCTrue().ParticleTree(), data);
}

template<typename T>
T& getHistogram(const data::ParticleTree_t& particletree,
                std::map<std::string, T>& perDecayHists,
                SmartHistFactory& HistFac) {
    const std::string& decaystring = utils::ParticleTools::GetDecayString(particletree);
    // search map only once even on insert
    auto it_h = perDecayHists.lower_bound(decaystring);
    if(it_h == perDecayHists.end() || perDecayHists.key_comp()(decaystring, it_h->first)) {
        // decaystring does not exist
        it_h = perDecayHists.emplace_hint(it_h, decaystring, T(HistFac, decaystring));
    }
    return it_h->second;
}

void EtapOmegaG::ProcessSig(const data::ParticleTree_t& particletree,
                            const data::Event::Data& data)
{
    auto& steps = sig_steps;

    const auto nParticles = data.Particles().GetAll().size();

    steps->Fill("Seen",1);

    if(nParticles != 5)
        return;

    steps->Fill("nParticles==5",1);

    const auto& photons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles().Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 4)
        return;
    steps->Fill("nPhotons==4",1);

    if(nProtons != 1)
        return;
    steps->Fill("nProtons==1",1);
    const data::ParticlePtr& proton = protons.front();

    if(!(proton->Candidate()->Detector() & Detector_t::Type_t::TAPS))
        return;
    steps->Fill("p in TAPS", 1);

    if(data.TriggerInfos().CBEenergySum() < 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    sig_perDecayHists_t& h = getHistogram(particletree, sig_perDecayHists, HistFac_sig);

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

    // proton coplanarity
    const double d_phi = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
    const interval<double> Proton_Copl_cut(-19, 19);
    if(!Proton_Copl_cut.Contains(d_phi))
        return;
    steps->Fill("Copl p in 2#sigma",1);
    h.Proton_Copl->Fill(d_phi);

    // bottom-up assignment of photons using Chi2

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

        const auto Chi2_Pi0 =  std_ext::sqr((tmp.Pi0.M() - Pi0.Mean) / Pi0.Sigma);
        const auto Chi2_Omega = std_ext::sqr((tmp.Omega.M() - Omega.Mean) / Omega.Sigma);
        const auto Chi2_EtaPrime = std_ext::sqr((tmp.EtaPrime.M() - EtaPrime_sig.Mean) / EtaPrime_sig.Sigma);

        tmp.Chi2 = Chi2_Pi0 + Chi2_Omega + Chi2_EtaPrime;

        h.Chi2_All->Fill(tmp.Chi2);

        if(tmp.Chi2<result.Chi2)
            result = tmp;
    }
    while(next_permutation(indices.begin(), indices.end(), comparer));

    // photon assignment successful?
    if(result.Chi2>10)
        return;
    steps->Fill("MinChi2<10",1);

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

void EtapOmegaG::ProcessRef(const data::ParticleTree_t& particletree,
                            const data::Event::Data& data)
{
    auto& steps = ref_steps;

    const auto nParticles = data.Particles().GetAll().size();

    steps->Fill("Seen",1);

    if(nParticles != 3)
        return;

    steps->Fill("nParticles==3",1);

    const auto& photons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles().Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 2)
        return;
    steps->Fill("nPhotons==2",1);

    if(nProtons != 1)
        return;
    steps->Fill("nProtons==1",1);
    const data::ParticlePtr& proton = protons.front();

    if(!(proton->Candidate()->Detector() & Detector_t::Type_t::TAPS))
        return;
    steps->Fill("p in TAPS", 1);


    if(data.TriggerInfos().CBEenergySum() < 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    ref_perDecayHists_t& h = getHistogram(particletree, ref_perDecayHists, HistFac_ref);

    // gamma combinatorics
    assert(photons.size() == 2);

    TLorentzVector photon_sum(0,0,0,0);
    for(const auto& p : photons) {
        photon_sum += *p;
    }
    const double gg_im = photon_sum.M();
    h.gg->Fill(gg_im);

    // proton coplanarity
    const double d_phi = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
    const interval<double> Proton_Copl_cut(-19, 19);
    if(!Proton_Copl_cut.Contains(d_phi))
        return;
    steps->Fill("Copl p in 2#sigma",1);
    h.Proton_Copl->Fill(d_phi);


    const double sigma = 2*EtaPrime_ref.Sigma;
    const interval<double> IM_etap_cut(EtaPrime_ref.Mean-sigma,EtaPrime_ref.Mean+sigma);

    if(!IM_etap_cut.Contains(gg_im))
        return;
    steps->Fill("IM #eta' in 2#sigma",1);

    h.IM_etap->Fill(photon_sum.M());

    // use tagged photon
    for(const auto& th : data.TaggerHits()) {
        /// \todo make prompt/random cut
        const TLorentzVector beam_target = th->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector mm = beam_target - photon_sum;
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
    canvas c_steps(GetName()+": Steps");
    c_steps << sig_steps << ref_steps << endc;

    for(const auto& it_map : sig_perDecayHists) {
        const sig_perDecayHists_t& h = it_map.second;
        if(h.IM_etap_omega->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it_map.first);
        c << h.gg << h.ggg << h.gggg << h.MM_gggg
          << h.Proton_Copl
          << h.Chi2_All << h.Chi2_Best
          << h.IM_pi0 << drawoption("colz") << h.IM_etap_omega
          << h.MM_etap
          << drawoption("colz") <<  h.Proton_ThetaPhi
          << h.Proton_Energy
          << endc;
    }

    for(const auto& it_map : ref_perDecayHists) {
        const ref_perDecayHists_t& h = it_map.second;
        if(h.IM_etap->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it_map.first);
        c << h.gg
          << h.Proton_Copl << h.IM_etap
          << h.MM_etap
          << drawoption("colz") << h.Proton_ThetaPhi
          << h.Proton_Energy
          << endc;
    }
}






AUTO_REGISTER_PHYSICS(EtapOmegaG)
