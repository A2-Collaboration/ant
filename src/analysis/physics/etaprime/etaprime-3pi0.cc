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



Etap3pi0::Etap3pi0(const std::string& name, PhysOptPtr opts) :
    Physics(name, opts),
      dataset(opts->GetOption("dataset"))
{
    BinSettings bs = BinSettings(1200);

    hNgammaMC  = HistFac.makeTH1D("# gamma MC true","# #gamma","# events",BinSettings(14));
    hNgamma    = HistFac.makeTH1D("# gamma","# #gamma","# events",BinSettings(14));
    h2g        = HistFac.makeTH1D("2 #gamma","2#gamma IM [MeV]","#",bs,"gg");
    h6g        = HistFac.makeTH1D("6 #gamma","6#gamma IM [MeV]","#",bs,"gggggg");

    ch_3pi0_IM_etap    = HistFac.makeTH1D("EtaPrime (3pi0)","EtaPrime IM [MeV]","events",bs,"ch_3pi0_IM_etap");
    ch_3pi0_IM_pi0     = HistFac.makeTH1D("Pi0 (3pi0)","Pi0 IM [MeV]","events",bs,"ch_3pi0_IM_pi0");

    ch_eta2pi0_IM_etap    = HistFac.makeTH1D("EtaPrime (eta2pi0)","EtaPrime IM [MeV]","events",bs,"ch_eta2pi0_IM_etap");
    ch_eta2pi0_IM_pions  = HistFac.makeTH1D("pions (eta2pi0)","Pi0 IM [MeV]","events",bs,"ch_eta2pi0_IM_pions");
    ch_eta2pi0_IM_etas  = HistFac.makeTH1D("etas (eta2pi0)","Pi0 IM [MeV]","events",bs,"ch_eta2pi0_IM_etas");

    dalitz_xy     = HistFac.makeTH2D("dalitz","x","y",BinSettings(100,0,0),BinSettings(100,0,0));
    dalitz_z      = HistFac.makeTH1D("dalitz - radial","z = x^{2} + y^{2}","#",BinSettings(100,0,0));

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

    result.Chi2_etaprime = std_ext::sqr( (result.etaprime.M() - 893.0) / 24.3);
    return result;
}

Etap3pi0::result_t Etap3pi0::MakeMC3pi0(const Event::Data& mcEvt )
{
    result_t result;

    result.Chi2 = 0;

    const auto& pions = mcEvt.Particles().Get(ParticleTypeDatabase::Pi0);
    const auto& etaprime = mcEvt.Particles().Get(ParticleTypeDatabase::EtaPrime);

    result.success = ( pions.size() == 3) && (etaprime.size() == 1 );

    if (result.success)
    {
        for ( unsigned i = 0 ; i < 3 ; ++i)
        {
            result.mesons[i].first = pions.at(i);
            result.mesons[i].second = 0;
        }
        result.etaprime = *(etaprime.at(0));
    }
    result.Chi2_etaprime = std_ext::sqr( (result.etaprime.M() - 893.0) / 24.3);

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

    result.Chi2_etaprime = std_ext::sqr( (result.etaprime.M() - 893.0) / 24.3);
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

    if (result_3pi0.Chi2_etaprime < 2)
    {
        result_3pi0.FillIm(ParticleTypeDatabase::Pi0, ch_3pi0_IM_pi0);
        result_3pi0.FillImEtaPrime(ch_3pi0_IM_etap);
    }

    if (result_eta2pi0.Chi2_etaprime < 2)
    {
        result_eta2pi0.FillIm(ParticleTypeDatabase::Pi0, ch_eta2pi0_IM_pions);
        result_eta2pi0.FillIm(ParticleTypeDatabase::Eta, ch_eta2pi0_IM_etas);
        result_eta2pi0.FillImEtaPrime(ch_eta2pi0_IM_etap);
    }

//    result_t result_mc = MakeMC3pi0(mcdata);
    DalitzVars channel(result_3pi0);
//    DalitzVars channel(result_mc);
//    DalitzVars refecence(result_eta2pi0);

    if (result_3pi0.Chi2_etaprime < 2 && result_3pi0.success)
    {
        dalitz_xy->Fill(channel.x,channel.y);
        dalitz_z->Fill(channel.z);
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
    canvas("Dalitz-Plots") << dalitz_z << drawoption("colz") << dalitz_xy << endc;
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



Etap3pi0::DalitzVars::DalitzVars(Etap3pi0::result_t r)
{
    s1 = (r.etaprime - *(r.mesons[0].first)).M2();
    s2 = (r.etaprime - *(r.mesons[1].first)).M2();
    s3 = (r.etaprime - *(r.mesons[2].first)).M2();

//    TMean = r.etaprime.M();
    for (const auto& meson: r.mesons)
        TMean -= meson.first->M();
    TMean = TMean / 3.0;

    x = ( r.mesons[0].first->Ek() - r.mesons[1].first->Ek() ) / ( TMath::Sqrt(3) * TMean);

    y = 0;
    for (const auto& meson: r.mesons)
        y += meson.first->M();
    y = y / ( 3 * r.etaprime.M());
    y = (y * r.mesons[2].first->Ek() / TMean) - 1;

    z = std_ext::sqr(x) + std_ext::sqr(y);

}

AUTO_REGISTER_PHYSICS(Etap3pi0)
