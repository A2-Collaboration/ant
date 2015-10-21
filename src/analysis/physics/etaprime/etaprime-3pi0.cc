#include "etaprime-3pi0.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "base/std_ext/math.h"
#include "data/Particle.h"

#include <algorithm>
#include <cassert>

#include "TTree.h"
#include "TCanvas.h"

using namespace std;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;



Etap3pi0::Etap3pi0(PhysOptPtr opts) :
    Physics("EtapOmegaG", opts),
      dataset(opts->GetOption("dataset"))
{
    BinSettings bs = BinSettings(1200);

    hNgammaMC  = HistFac.makeTH1D("# gamma MC true","# #gamma","# events",BinSettings(14));
    hNgamma    = HistFac.makeTH1D("# gamma","# #gamma","# events",BinSettings(14));
    h2g        = HistFac.makeTH1D("2 #gamma","2#gamma IM [MeV]","#",bs,"gg");
    h6g        = HistFac.makeTH1D("6 #gamma","6#gamma IM [MeV]","#",bs,"gggggg");

    ch_3pi0_IM_etap    = HistFac.makeTH1D("EtaPrime (3pi0)","EtaPrime IM [MeV]","events",bs,"IM_etap");
    ch_3pi0_IM_pi0     = HistFac.makeTH1D("Pi0 (3pi0)","Pi0 IM [MeV]","events",bs,"IM_pi0");

    ch_eta2pi0_IM_etap    = HistFac.makeTH1D("EtaPrime (eta2pi0)","EtaPrime IM [MeV]","events",bs,"IM_etap");
    ch_eta2pi0_IM_pions  = HistFac.makeTH1D("pions (eta2pi0)","Pi0 IM [MeV]","events",bs,"IM_pi0");
    ch_eta2pi0_IM_etas  = HistFac.makeTH1D("etas (eta2pi0)","Pi0 IM [MeV]","events",bs,"IM_pi0");

    //dalitz     = HistFac.makeTH2D("dalitz","IM^{2}(#pi^{0} _{1} + #pi^{0} _{2}) [GeV^{2}]","IM^{2}(#pi^{0} _{1} + #pi^{0} _{3}) [GeV^{2}]",BinSettings(100,0,1000),BinSettings(100,0,1000));

    tree       = new TTree("data","data");


    pi01.SetBranches(tree,"pi01");
    pi02.SetBranches(tree,"pi02");
    pi03.SetBranches(tree,"pi03");
    MMproton.SetBranches(tree,"MMproton");
    etaprime.SetBranches(tree,"etaprime");

    tree->Branch("imsqr12",&imsqrP12);
    tree->Branch("imsqr13",&imsqrP13);
    tree->Branch("imsqr23",&imsqrP23);

}

void Etap3pi0::FillCrossChecks(const ParticleList& photons, const ParticleList& mcphotons)
{
    hNgamma->Fill(photons.size());
    hNgammaMC->Fill(mcphotons.size());

    if (photons.size() != 6)
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
}

Etap3pi0::result_t Etap3pi0::Make3pi0(const ParticleList& photons)
{
    result_t result;

    for ( const auto& pairs: combinations)
    {

        result_t tmp;
        tmp.Chi2 = 0;

        for(unsigned i=0;i<pairs.size();i++)
        {
            tmp.g_final[pairs.at(i).first] = photons[2*i];
            tmp.g_final[pairs.at(i).second] = photons[2*i+1];
        }

        for(unsigned i=0;i<tmp.mesons.size();i++)
        {
            tmp.mesons[i].first = make_shared<Particle>(ParticleTypeDatabase::Pi0, *(tmp.g_final[2*i]) + *(tmp.g_final[2*i+1]));
            tmp.Chi2 += std_ext::sqr((tmp.mesons[i].first->M() - 126) / 15); // width and center from fit
        }
        if(tmp.Chi2<result.Chi2)
            result = move(tmp);
    }

    result.success = true;
    for (const auto& p: result.mesons )
    {
//        ch_3pi0_IM_pi0->Fill(p->M());
        result.etaprime += *(p.first);
    }
//    ch_3pi0_IM_etap->Fill(result.etaprime.M());

    return result;
}

Etap3pi0::result_t Etap3pi0::MakeEta2pi0(const ParticleList& photons)
{
    result_t result;

    for ( const auto& pairs: combinations)
    {

        result_t tmp;
        tmp.Chi2 = 0;

        for(unsigned i=0;i<pairs.size();i++)
        {
            tmp.g_final[pairs.at(i).first] = photons[2*i];
            tmp.g_final[pairs.at(i).second] = photons[2*i+1];
        }


        for (unsigned etaIndex = 0 ; etaIndex < tmp.mesons.size() ; ++etaIndex)
        {
            tmp.mesons[etaIndex].first = make_shared<Particle>(ParticleTypeDatabase::Eta,*(tmp.g_final[2*etaIndex]) + *(tmp.g_final[2*etaIndex+1]));
            tmp.Chi2 =  std_ext::sqr((tmp.mesons[etaIndex].first->M() - 515.5) / 19.4);        // width and center from fit

            unsigned piIndex = ( etaIndex + 1 ) % 3;
            tmp.mesons[piIndex].first = make_shared<Particle>(ParticleTypeDatabase::Pi0,*(tmp.g_final[2*piIndex]) + *(tmp.g_final[2*piIndex+1]));
            tmp.Chi2 += std_ext::sqr((tmp.mesons[piIndex].first->M() - 126) / 15);

            piIndex = ( etaIndex + 2 ) % 3;
            tmp.mesons[piIndex].first = make_shared<Particle>(ParticleTypeDatabase::Pi0,*(tmp.g_final[2*piIndex]) + *(tmp.g_final[2*piIndex+1]));
            tmp.Chi2 += std_ext::sqr((tmp.mesons[piIndex].first->M() - 126) / 15);

            if(tmp.Chi2<result.Chi2)
                result = move(tmp);
        }

    }

    result.success = true;
    for (const auto& p: result.mesons )
        result.etaprime += *(p.first);

    /*
    ch_eta2pi0_IM_etas->Fill(result.mesons[etafound]->M());
    ch_eta2pi0_IM_pions->Fill(result.mesons[( etafound + 1 ) % 3]->M());
    ch_eta2pi0_IM_pions->Fill(result.mesons[( etafound + 2 ) % 3]->M());
    ch_eta2pi0_IM_etap->Fill(result.etaprime.M());
    */

    return result;
}


void Etap3pi0::ProcessEvent(const data::Event& event)
{
    const auto& data   = event.Reconstructed();
    const auto& mcdata = event.MCTrue();

    const auto& photons   = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& mcphotons = mcdata.Particles().Get(ParticleTypeDatabase::Photon);

    FillCrossChecks(photons,mcphotons);

    if (photons.size() != 6 )
        return;

    result_t result_3pi0    = Make3pi0(photons);
    result_t result_eta2pi0 = MakeEta2pi0(photons);

    pi01     = ParticleVars(*(result_3pi0.mesons[0].first));
    pi02     = ParticleVars(*(result_3pi0.mesons[1].first));
    pi03     = ParticleVars(*(result_3pi0.mesons[2].first));
    etaprime = ParticleVars(result_3pi0.etaprime, ParticleTypeDatabase::EtaPrime);

    if (result_3pi0.Chi2 < 2)
    {
        result_3pi0.FillIm(ParticleTypeDatabase::Pi0, ch_3pi0_IM_pi0);
        result_3pi0.FillImEtaPrime(ch_3pi0_IM_etap);
    }

    if (result_eta2pi0.Chi2 < 2)
    {
        result_eta2pi0.FillIm(ParticleTypeDatabase::Pi0, ch_eta2pi0_IM_pions);
        result_eta2pi0.FillIm(ParticleTypeDatabase::Eta, ch_eta2pi0_IM_etas);
        result_eta2pi0.FillImEtaPrime(ch_eta2pi0_IM_etap);
    }



  //TLorentzVector sum;
  //sum = result.Pi0[0] + result.Pi0[1];
  //imsqrP12 = sum.M2();

  //sum = result.Pi0[0] + result.Pi0[2];
  //imsqrP13 = sum.M2();

  //sum = result.Pi0[1] + result.Pi0[2];
  //imsqrP23 = sum.M2();

  // for(const auto& taggerhit : data.TaggerHits())
  // {
  //     const TLorentzVector beam_target = taggerhit->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
  //     const Particle missing(ParticleTypeDatabase::Proton, beam_target - etap);
  //     MMproton = ParticleVars(missing);
  // }

    tree->Fill();

    //dalitz->Fill(imsqrP12/1000.0,imsqrP13/1000.0);
}


void Etap3pi0::ShowResult()
{
    canvas("Crosschecks")       << h2g << h6g
                                << hNgamma << hNgammaMC
                                << endc;
    canvas("Invaraiant Masses") << ch_eta2pi0_IM_pions << ch_eta2pi0_IM_etas << ch_eta2pi0_IM_etap
                                << ch_3pi0_IM_pi0 << ch_3pi0_IM_etap
                                //<< IM_proton << IM_mmproton
                                << endc;
    /*
    canvas("Dalitz-Plots") << drawoption("colz") << dalitz << endc;
    */
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



void Etap3pi0::result_t::FillIm(const ant::ParticleTypeDatabase::Type& type, TH1D* hist)
{
    for (const auto& meson: mesons)
    {
        if ( meson.first->Type() == type )
            hist->Fill(meson.first->M());
    }
}

void Etap3pi0::result_t::FillImEtaPrime(TH1D* hist)
{
    hist->Fill(etaprime.M());
}

AUTO_REGISTER_PHYSICS(Etap3pi0, "Etap3pi0")
