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




EtapOmegaG::EtapOmegaG(const std::string& name, PhysOptPtr opts) : Physics(name, opts),
    HistFac_sig("Sig",HistFac),
    HistFac_ref("Ref",HistFac),
    sig_hists(HistFac_sig),
    ref_hists(HistFac_ref)
{

    treeSig = HistFac_sig.makeTTree("treeSig");
    treeRef = HistFac_ref.makeTTree("treeRef");

    sig_perDecayHists.emplace_back(
                "Signal",
                HistFac_sig,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_2pi0",
                HistFac_sig,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct2Pi0_4g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_3pi0",
                HistFac_sig,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct3Pi0_6g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_OmegaPi0g",
                HistFac_sig,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_1pi0",
                HistFac_sig,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct1Pi0_2g)
                );
    sig_perDecayHists.emplace_back("Signal_Bkg_Other", HistFac_sig);


    ref_perDecayHists.emplace_back(
                "Reference",
                HistFac_ref,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g)
                );
    ref_perDecayHists.emplace_back(
                "Ref_Bkg_2pi0",
                HistFac_ref,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct2Pi0_4g)
                );
    ref_perDecayHists.emplace_back(
                "Ref_Bkg_1pi0",
                HistFac_ref,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct1Pi0_2g)
                );
    ref_perDecayHists.emplace_back("Ref_Bkg_Other", HistFac_ref);
}

EtapOmegaG::histogram_t::histogram_t(SmartHistFactory histFac)
{
    Steps = histFac.makeTH1D("Steps", "", "#", BinSettings(10),"steps");
    MissedBkg = histFac.makeTH1D("Missed Bkg channels", "", "#", BinSettings(20),"missed_bkg");
}

EtapOmegaG::sig_perDecayHists_t::sig_perDecayHists_t(SmartHistFactory HistFac)
{
    BinSettings bins_im(1200);

    gggg = HistFac.makeTH1D("4#gamma IM","4#gamma IM / MeV","events",bins_im,"gggg");
    ggg = HistFac.makeTH1D("3#gamma IM","3#gamma IM / MeV","events",bins_im,"ggg");
    gg = HistFac.makeTH1D("2#gamma IM","2#gamma IM / MeV","events",bins_im,"gg");

    IM_gg_gg = HistFac.makeTH2D("IM 2#gamma vs. 2#gamma (Goldhaber plot)",
                                "2#gamma IM / MeV", "2#gamma IM / MeV",
                                BinSettings(400, 0, 800),
                                BinSettings(400, 0, 800),
                                "IM_gg_gg"
                                );

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
    Chi2_Best = HistFac.makeTH1D("#chi^{2} minimal","#chi^{2}","",BinSettings(300,0,12),"Chi2_Min");

    Proton_ThetaPhi = HistFac.makeTH2D("p #delta(#theta-#phi)",
                                       "#delta#theta / degree",
                                       "#delta#phi / degree",
                                       BinSettings(100, -10, 10),
                                       BinSettings(100, -30, 30),
                                       "Proton_ThetaPhi"
                                       );
    Proton_Energy = HistFac.makeTH1D("p #delta(E)","#deltaE / MeV","",BinSettings(400,-50,350),"Proton_Energy");
}

EtapOmegaG::ref_perDecayHists_t::ref_perDecayHists_t(SmartHistFactory HistFac)
{
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
const T& getHistogram(const data::ParticleTree_t& particletree,
                      const std::vector<EtapOmegaG::perDecayHists_t<T>>& perDecayHists,
                      bool& other
                      ) {
    assert(!perDecayHists.empty());
    other = true;
    if(!particletree)
        return perDecayHists.back().PerDecayHists;
    for(size_t i=0;i<perDecayHists.size()-1;i++) {
        auto& item = perDecayHists[i];
        if(particletree->IsEqual(item.Tree, utils::ParticleTools::MatchByParticleName)) {
            other = false;
            return item.PerDecayHists;
        }
    }
    return perDecayHists.back().PerDecayHists;
}

void EtapOmegaG::ProcessSig(const data::ParticleTree_t& particletree,
                            const data::Event::Data& data)
{
    auto& steps = sig_hists.Steps;

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

    bool other_channel = false;
    const sig_perDecayHists_t& h = getHistogram(particletree, sig_perDecayHists, other_channel);

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

    // fill Goldhaber plot
    for(auto i : std::initializer_list<std::vector<unsigned>>({{0,1,2,3},{0,2,1,3},{0,3,1,2}})) {
        auto& p = photons;
        h.IM_gg_gg->Fill((*p[i[0]]+*p[i[1]]).M(), (*p[i[2]]+*p[i[3]]).M());
    }

    //h.IM_gg_gg->Fill(  );

    // was this some unidentified channel?
    if(other_channel) {
        sig_hists.MissedBkg->Fill(utils::ParticleTools::GetDecayString(particletree).c_str(),1);
    }

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
    auto& steps = ref_hists.Steps;

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

    bool other_channel = false;
    const ref_perDecayHists_t& h = getHistogram(particletree, ref_perDecayHists, other_channel);

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

    // was this some unidentified channel?
    if(other_channel) {
        ref_hists.MissedBkg->Fill(utils::ParticleTools::GetDecayString(particletree).c_str(),1);
    }

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
}

void EtapOmegaG::ShowResult()
{
    canvas c_steps(GetName()+": Steps");
    c_steps << sig_hists.Steps << sig_hists.MissedBkg
            << ref_hists.Steps << ref_hists.MissedBkg
            << endc;

    for(const auto& it : sig_perDecayHists) {
        const sig_perDecayHists_t& h = it.PerDecayHists;
        if(h.IM_etap_omega->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it.ShortName);
        c << h.gg << h.ggg << h.gggg << h.MM_gggg
          << h.Proton_Copl
          << h.Chi2_All << h.Chi2_Best
          << h.IM_pi0
          << drawoption("colz") << h.IM_etap_omega
          << drawoption("colz") << h.IM_gg_gg
          << h.MM_etap
          << drawoption("colz") <<  h.Proton_ThetaPhi
          << h.Proton_Energy
          << endc;
    }

    for(const auto& it : ref_perDecayHists) {
        const ref_perDecayHists_t& h = it.PerDecayHists;
        if(h.IM_etap->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it.ShortName);
        c << h.gg
          << h.Proton_Copl << h.IM_etap
          << h.MM_etap
          << drawoption("colz") << h.Proton_ThetaPhi
          << h.Proton_Energy
          << endc;
    }
}










AUTO_REGISTER_PHYSICS(EtapOmegaG)
