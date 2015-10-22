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

    mcdalitz_xy     = HistFac.makeTH2D("dalitz mc","s1","s3",BinSettings(100,0,0),BinSettings(100,0,0));
    mcdalitz_z      = HistFac.makeTH1D("dalitz - radial mc","z = x^{2} + y^{2}","#",BinSettings(100,0,0));

    dalitz_xy     = HistFac.makeTH2D("dalitz","s1","s3",BinSettings(100,0,0),BinSettings(100,0,0));
    dalitz_z      = HistFac.makeTH1D("dalitz - radial","z = x^{2} + y^{2}","#",BinSettings(100,0,0));

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
        tmp.Chi2_intermediate = 0;

        for(unsigned i=0;i<pairs.size();i++)
        {
            tmp.g_final[pairs.at(i).first] = photons[2*i];
            tmp.g_final[pairs.at(i).second] = photons[2*i+1];
        }

        for(unsigned i=0;i<tmp.mesons.size();i++)
        {
            tmp.mesons[i].first = make_shared<Particle>(ParticleTypeDatabase::Pi0, *(tmp.g_final[2*i]) + *(tmp.g_final[2*i+1]));
            tmp.Chi2_intermediate += std_ext::sqr((tmp.mesons[i].first->M() - 126) / 15); // width and center from fit
        }
        if(tmp.Chi2_intermediate<result.Chi2_intermediate)
            result = move(tmp);
    }

    result.success = true;
    for (const auto& p: result.mesons )
    {
        result.mother += *(p.first);
    }

    result.Chi2_mother = std_ext::sqr( (result.mother.M() - 893.0) / 24.3);
    return result;
}

Etap3pi0::result_t Etap3pi0::MakeMC3pi0(const Event::Data& mcEvt )
{
    result_t result;

    result.Chi2_intermediate = 0;
    result.Chi2_mother = 0;

    const auto& pions = mcEvt.Intermediates().Get(ParticleTypeDatabase::Pi0);
    const auto& etaprime = mcEvt.Intermediates().Get(ParticleTypeDatabase::EtaPrime);

    result.success = ( pions.size() == 3) && (etaprime.size() == 1 );

    if (result.success)
    {
        for ( unsigned i = 0 ; i < 3 ; ++i)
        {
            result.mesons[i].first = pions.at(i);
            result.mesons[i].second = 0;
        }
        result.mother = *(etaprime.at(0));
    }

    return result;
}

Etap3pi0::result_t Etap3pi0::MakeEta2pi0(const ParticleList& photons)
{
    result_t result;

    for ( const auto& pairs: combinations)
    {

        result_t tmp;
        tmp.Chi2_intermediate = 0;

        for(unsigned i=0;i<pairs.size();i++)
        {
            tmp.g_final[pairs.at(i).first] = photons[2*i];
            tmp.g_final[pairs.at(i).second] = photons[2*i+1];
        }


        for (unsigned etaIndex = 0 ; etaIndex < tmp.mesons.size() ; ++etaIndex)
        {
            tmp.mesons[etaIndex].first = make_shared<Particle>(ParticleTypeDatabase::Eta,*(tmp.g_final[2*etaIndex]) + *(tmp.g_final[2*etaIndex+1]));
            tmp.Chi2_intermediate =  std_ext::sqr((tmp.mesons[etaIndex].first->M() - 515.5) / 19.4);        // width and center from fit

            unsigned piIndex = ( etaIndex + 1 ) % 3;
            tmp.mesons[piIndex].first = make_shared<Particle>(ParticleTypeDatabase::Pi0,*(tmp.g_final[2*piIndex]) + *(tmp.g_final[2*piIndex+1]));
            tmp.Chi2_intermediate += std_ext::sqr((tmp.mesons[piIndex].first->M() - 126) / 15);

            piIndex = ( etaIndex + 2 ) % 3;
            tmp.mesons[piIndex].first = make_shared<Particle>(ParticleTypeDatabase::Pi0,*(tmp.g_final[2*piIndex]) + *(tmp.g_final[2*piIndex+1]));
            tmp.Chi2_intermediate += std_ext::sqr((tmp.mesons[piIndex].first->M() - 126) / 15);

            if(tmp.Chi2_intermediate<result.Chi2_intermediate)
                result = move(tmp);
        }

    }

    result.success = true;
    for (const auto& p: result.mesons )
        result.mother += *(p.first);

    result.Chi2_mother = std_ext::sqr( (result.mother.M() - 906.0) / 26.3);
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
    result_t result_mc      = MakeMC3pi0(mcdata);


    const double chi2cut(3);

    if (result_3pi0.Chi2_mother < chi2cut && result_3pi0.Chi2_intermediate < chi2cut)
    {
        FillIm(result_3pi0, ParticleTypeDatabase::Pi0, ch_3pi0_IM_pi0);
        FillImEtaPrime(result_3pi0,ch_3pi0_IM_etap);
    }

    if (result_eta2pi0.Chi2_mother < chi2cut && result_eta2pi0.Chi2_intermediate < chi2cut)
    {
        FillIm(result_eta2pi0, ParticleTypeDatabase::Pi0, ch_eta2pi0_IM_pions);
        FillIm(result_eta2pi0, ParticleTypeDatabase::Eta, ch_eta2pi0_IM_etas);
        FillImEtaPrime(result_eta2pi0, ch_eta2pi0_IM_etap);
    }

    if ( result_mc.success)
    {
        DalitzVars channel(result_mc);
        mcdalitz_xy->Fill(channel.s1,channel.s3);
        mcdalitz_z->Fill(channel.z);
    }
    if (result_3pi0.Chi2_mother < chi2cut && result_3pi0.Chi2_intermediate < chi2cut )
    {
        DalitzVars channel(result_3pi0);
        dalitz_xy->Fill(channel.s1,channel.s3);
        dalitz_z->Fill(channel.z);
    }

}


void Etap3pi0::ShowResult()
{
    canvas("Crosschecks")       << h2g << h6g
                                << hNgamma << hNgammaMC
                                << endc;

    canvas("Invaraiant Masses") << ch_eta2pi0_IM_pions << ch_eta2pi0_IM_etas << ch_eta2pi0_IM_etap
                                << ch_3pi0_IM_pi0 << ch_3pi0_IM_etap
                                << endc;

    canvas("Dalitz-Plots")      << drawoption("colz") << dalitz_xy << mcdalitz_xy << endc;
}



void Etap3pi0::FillIm(const Etap3pi0::result_t& result, const ant::ParticleTypeDatabase::Type& type, TH1D* hist)
{
    for (const auto& meson: result.mesons)
    {
        if ( meson.first->Type() == type )
            hist->Fill(meson.first->M());
    }
}

void Etap3pi0::FillImEtaPrime(const Etap3pi0::result_t& result, TH1D* hist)
{
    hist->Fill(result.mother.M());
}



Etap3pi0::DalitzVars::DalitzVars(Etap3pi0::result_t r)
{
    s1 = (r.mother - *(r.mesons[0].first)).M2();
    s2 = (r.mother - *(r.mesons[1].first)).M2();
    s3 = (r.mother - *(r.mesons[2].first)).M2();

//    TMean = r.etaprime.M();
    for (const auto& meson: r.mesons)
        TMean -= meson.first->M();
    TMean = TMean / 3.0;

    x = ( r.mesons[0].first->Ek() - r.mesons[1].first->Ek() ) / ( TMath::Sqrt(3) * TMean);

    y = 0;
    for (const auto& meson: r.mesons)
        y += meson.first->M();
    y = y / ( 3 * r.mother.M());
    y = (y * r.mesons[2].first->Ek() / TMean) - 1;

    z = std_ext::sqr(x) + std_ext::sqr(y);

}

AUTO_REGISTER_PHYSICS(Etap3pi0)
