#include "etaprime-3pi0.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "base/std_ext/math.h"
#include "data/Particle.h"
#include "utils/particle_tools.h"

#include <algorithm>
#include <cassert>
#include <chrono>

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
    BinSettings bs_im = BinSettings(1200);
    string cat("hist");

    // no category
    cat = "proton";
    AddHist1D(cat,"ProtonCandidateAngles",          "Proton Candidate Angles", "#Theta [#circ]", "#", BinSettings(180));
    AddHist1D(cat,"mcProtonAngles",                 "MC Proton Angles", "#Theta [#circ]", "#", BinSettings(180));


    cat = "P";
    AddHist2D(cat, "all", "All", "#Theta [#circ]","E [MeV]",BinSettings(360,0,180),bs_im);
    AddHist2D(cat, "gamma_all", "All #gamma", "#Theta [#circ]","E [MeV]",BinSettings(360,0,180),bs_im);
    AddHist2D(cat, "gamma_6", "All 6 #gamma", "#Theta [#circ]","E [MeV]",BinSettings(360,0,180),bs_im);
    AddHist2D(cat, "gamma_signal", "signal 6 #gamma", "#Theta [#circ]","E [MeV]",BinSettings(360,0,180),bs_im);
    AddHist2D(cat, "gamma_ref", "reference 6 #gamma", "#Theta [#circ]","E [MeV]",BinSettings(360,0,180),bs_im);



    cat = "xc";             // crosschecks
    AddHist1D(cat, "NgammaMC", "# gamma MC true","# #gamma","# events",BinSettings(14));
    AddHist1D(cat, "Ngamma"  , "# gamma","# #gamma","# events",BinSettings(14));
    AddHist1D(cat, "IM_2g"   , "2 #gamma","2#gamma IM [MeV]","#",bs_im);
    AddHist1D(cat, "IM_6g"   , "6 #gamma","6#gamma IM [MeV]","#",bs_im);

    AddHist1D(cat, "NTagger" , "# Taggerhits","#","",BinSettings(3));
    AddHist1D(cat, "NProtons", "# Protons",   "#","",BinSettings(6));

    cat = "chi2s";
    AddHist1D(cat,"signal_pi0","#chi^{2} of #pi{0} candidates (signal)","chi^{2}","#",BinSettings(500,0,50));
    AddHist1D(cat,"signal_intermediate","#chi^{2} for selection (signal)","chi^{2}","#",BinSettings(100,0,5));
    AddHist1D(cat,"signal_etaprime","#chi^{2} for #eta' (signal)","chi^{2}","#",BinSettings(100,0,5));

    AddHist1D(cat,"ref_pi0","#chi^{2} of #pi^{0} candidates (reference)","chi^{2}","#",BinSettings(500,0,50));
    AddHist1D(cat,"ref_eta","#chi^{2} of #eta candidates (reference)","chi^{2}","#",BinSettings(500,0,50));
    AddHist1D(cat,"ref_intermediate","#chi^{2} for selection (reference)","chi^{2}","#",BinSettings(100,0,5));
    AddHist1D(cat,"ref_etaprime","#chi^{2} for #eta' (reference)","chi^{2}","#",BinSettings(100,0,5));

    /*
    cat = "kinfit";
    AddHist1D(cat,"signal_niter","","# iterations","#",BinSettings(100,0,5));
    AddHist1D(cat,"signal_chi2","#chi^{2} for kinfit (signal)","chi^{2}","#",BinSettings(100,0,5));
    AddHist1D(cat,"signal_egamma_before", "Photon Energy before (signal)", "E_{#gamma} before [MeV] (signal)","#",BinSettings(100,0,5));
    AddHist1D(cat,"signal_egamma_after", "Photon Energy after (signal)", "E_{#gamma} after [MeV] ","#",BinSettings(100,0,5));
    */


    cat = "signal";         //signal: eta' --> 3 pi0
    AddHist1D(cat,"IM_etap"    , "EtaPrime (3pi0)","EtaPrime IM [MeV]","events",bs_im);
    AddHist1D(cat,"IM_pi0"     , "Pi0 (3pi0)","Pi0 IM [MeV]","events",bs_im);


    cat = "ref";            //reference: eta' --> eta 2 pi0
    AddHist1D(cat,"IM_etap"  , "EtaPrime (eta2pi0)","EtaPrime IM [MeV]","events",bs_im);
    AddHist1D(cat,"IM_pions" , "pions (eta2pi0)","Pi0 IM [MeV]","events",bs_im);
    AddHist1D(cat,"IM_etas"  , "etas (eta2pi0)","Pi0 IM [MeV]","events",bs_im);

    cat = "dalitz";
    AddHist2D(cat, "mc_xy", "dalitz mc","s1","s3",BinSettings(100,0,0),BinSettings(100,0,0));
    AddHist1D(cat, "mc_z",  "dalitz - radial mc","z = x^{2} + y^{2}","#",BinSettings(100,0,0));
    AddHist2D(cat, "xy",    "dalitz","s1","s3",BinSettings(100,0,0),BinSettings(100,0,0));
    AddHist1D(cat, "z",     "dalitz - radial","z = x^{2} + y^{2}","#",BinSettings(100,0,0));

    cat = "steps";
    AddHist1D(cat, "evcount", "events after steps", "", "# events", BinSettings(5));

    cat = "channels";
    AddHist1D(cat,"nocut",              "6 #gamma, no cut", "", "#", BinSettings(15));
    AddHist1D(cat,"signal_chi2",        "6 #gamma, #chi^{2} cut (signal)", "", "#", BinSettings(15));
    AddHist1D(cat,"ref_chi2",           "6 #gamma, #chi^{2} cut (reference)", "", "#", BinSettings(15));
    AddHist1D(cat,"mc_true",            "mc true for signal, ref, bkg", "", "#", BinSettings(3));

    signal_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_3Pi0_6g);
    reference_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g);
    bkg_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct3Pi0_6g);
}

void Etap3pi0::FillCrossChecks(const ParticleList& photons, const ParticleList& mcphotons)
{
    hists.at("xc").at("Ngamma")->Fill(photons.size());
    hists.at("xc").at("NgammaMC")->Fill(mcphotons.size());

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


    fill_combinations(hists.at("xc").at("IM_2g"), 2, photons);
    fill_combinations(hists.at("xc").at("IM_6g"), 6, photons);
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

Etap3pi0::result_t Etap3pi0::Make3pi0(const ParticleList& photons)
{
    result_t result;

    const double pi0Chi2Cut = 3;

    bool found = false;
    for ( const auto& pairs: combinations)
    {

        result_t tmp;
        tmp.Chi2_intermediate = 0;

        for(unsigned i=0;i<pairs.size();i++)
        {
            tmp.g_final[pairs.at(i).first] = photons[2*i];
            tmp.g_final[pairs.at(i).second] = photons[2*i+1];
        }

        bool skip_comination = false;
        for(unsigned i=0;i<tmp.mesons.size();i++)
        {
            auto pi0candidate = make_shared<Particle>(ParticleTypeDatabase::Pi0, *(tmp.g_final[2*i]) + *(tmp.g_final[2*i+1]));
            double chi2 = std_ext::sqr((pi0candidate->M() - 126) / 15); // width and center from fit
            hists.at("chi2s").at("signal_pi0")->Fill(chi2);
            if ( chi2 > pi0Chi2Cut )
            {
                skip_comination = true;
                break;
            }
            tmp.mesons[i].first = pi0candidate;
            tmp.Chi2_intermediate += chi2;
        }
        if (!skip_comination)
        {
            if(tmp.Chi2_intermediate<result.Chi2_intermediate)
            {
                result = move(tmp);
                found = true;
            }
        }
    }

    result.success = found;
    if (found)
    {
        hists.at("chi2s").at("signal_intermediate")->Fill(result.Chi2_intermediate);
        for (const auto& p: result.mesons )
        {
            result.mother += *(p.first);
        }

        result.Chi2_mother = std_ext::sqr( (result.mother.M() - 893.0) / 24.3);
        hists.at("chi2s").at("signal_etaprime")->Fill(result.Chi2_mother);
    }
    return result;
}
Etap3pi0::result_t Etap3pi0::MakeEta2pi0(const ParticleList& photons)
{
    result_t result;
    const double pi0Chi2cut = 3;
    const double etaChi2cut = 3;
    bool found = false;

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
            unsigned etacount = 0;
            unsigned picount = 0;

            auto etacandidate = make_shared<Particle>(ParticleTypeDatabase::Eta,*(tmp.g_final[2*etaIndex]) + *(tmp.g_final[2*etaIndex+1]));
            double chi2eta    = std_ext::sqr((etacandidate->M() - 515.5) / 19.4);        // width and center from fit
            hists.at("chi2s").at("ref_eta")->Fill(chi2eta);
            if ( chi2eta < etaChi2cut)
                etacount++;

            unsigned pi1Index = ( etaIndex + 1 ) % 3;
            auto pi1candidate = make_shared<Particle>(ParticleTypeDatabase::Pi0,*(tmp.g_final[2*pi1Index]) + *(tmp.g_final[2*pi1Index+1]));
            double chipi1     = std_ext::sqr((pi1candidate->M() - 126) / 15);
            hists.at("chi2s").at("ref_pi0")->Fill(chipi1);
            if ( chipi1 < pi0Chi2cut )
                picount++;

            unsigned pi2Index = ( etaIndex + 2 ) % 3;
            auto pi2candidate =  make_shared<Particle>(ParticleTypeDatabase::Pi0,*(tmp.g_final[2*pi2Index]) + *(tmp.g_final[2*pi2Index+1]));
            double chipi2     =  std_ext::sqr((pi2candidate->M() - 126) / 15);
            hists.at("chi2s").at("ref_pi0")->Fill(chipi2);
            if (chipi2 < pi0Chi2cut)
                picount++;

            if ( picount == 2 && etacount == 1 )
            {
                found = true;
                tmp.mesons[etaIndex].first = etacandidate;
                tmp.Chi2_intermediate =  chi2eta;
                tmp.mesons[pi1Index].first = pi1candidate;
                tmp.Chi2_intermediate += chipi1;
                tmp.mesons[pi2Index].first = pi2candidate;
                tmp.Chi2_intermediate += chipi2;

                if(tmp.Chi2_intermediate<result.Chi2_intermediate)
                    result = move(tmp);
            }
        }

    }


    result.success  = found;
    if (found)
    {
        hists.at("chi2s").at("ref_intermediate")->Fill(result.Chi2_intermediate);
        for (const auto& p: result.mesons )
            result.mother += *(p.first);

        result.Chi2_mother = std_ext::sqr( (result.mother.M() - 906.0) / 26.3);
        hists.at("chi2s").at("ref_etaprime")->Fill(result.Chi2_mother);
    }
    return result;
}

bool Etap3pi0::MakeMCProton(const data::Event::Data& mcdata, ParticlePtr& proton)
{
   const auto& protonlist = mcdata.Particles().Get(ParticleTypeDatabase::Proton);
   if (protonlist.size() != 1)
       return false;
   proton = protonlist.at(0);
   return true;
}

void Etap3pi0::ProcessEvent(const data::Event& event)
{
    const auto& data   = event.Reconstructed();
    const auto& mcdata = event.MCTrue();

    if ( mcdata.ParticleTree() )
    {
        if (mcdata.ParticleTree()->IsEqual(signal_tree, utils::ParticleTools::MatchByParticleName))
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree()).c_str(),1);
        if (mcdata.ParticleTree()->IsEqual(reference_tree, utils::ParticleTools::MatchByParticleName))
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree()).c_str(),1);
        if (mcdata.ParticleTree()->IsEqual(bkg_tree, utils::ParticleTools::MatchByParticleName))
            hists.at("channels").at("mc_true")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree()).c_str(),1);
    }

    for (const auto& p: data.Particles().GetAll())
        hists["P"]["all"]->Fill(p->Theta()*TMath::RadToDeg(),p->Ek());

    const auto& photons            = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& protonCandidates   = data.Particles().Get(ParticleTypeDatabase::Proton);
    const auto& mcprotons          = mcdata.Particles().Get(ParticleTypeDatabase::Proton);

    for (const auto& ph: photons)
        hists["P"]["gamma_all"]->Fill(ph->Theta()*TMath::RadToDeg(),ph->Ek());

    hists.at("steps").at("evcount")->Fill("all",1);

    // MC Proton
    if (mcprotons.size() != 1)
        return;
    hists.at("steps").at("evcount")->Fill("req. mc proton",1);

    auto mcproton = mcprotons.at(0);
    hists.at("proton").at("mcProtonAngles")->Fill(mcproton->Theta() * TMath::RadToDeg());
    if (mcproton->Theta() * TMath::RadToDeg() > 20 )
        return;
    hists.at("steps").at("evcount")->Fill("mc proton angle < 20",1);



    const auto& mcphotons = mcdata.Particles().Get(ParticleTypeDatabase::Photon);

    FillCrossChecks(photons,mcphotons);
    hists.at("xc").at("NTagger")->Fill(data.TaggerHits().size());

    if (photons.size() != 6 )
        return;

    for (const auto& ph: photons)
        hists["P"]["gamma_6"]->Fill(ph->Theta()*TMath::RadToDeg(),ph->Ek());

    hists.at("steps").at("evcount")->Fill("req. 6 #gamma",1);
    hists.at("channels").at("nocut")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree()).c_str(),1);

    // we need a tagger hit
    if (data.TaggerHits().size() != 1)
        return;
    hists.at("steps").at("evcount")->Fill("req. 1 tagger hit",1);
    vector<ParticlePtr> protons;
    // only take proton if in TAPS and only select single
    for ( const auto& pcandidate: protonCandidates)
    {
        double thetaAngle = pcandidate->Theta() * TMath::RadToDeg() ;
        hists.at("proton").at("ProtonCandidateAngles")->Fill(thetaAngle);
        if ( thetaAngle < 20)
        {
            protons.push_back(pcandidate);
            continue;
        }
        return;
    }
    hists.at("steps").at("evcount")->Fill("proton angle < 20",1);
    hists.at("xc").at("NProtons")->Fill(protons.size());
    if (protons.size() == 0)
        hists.at("steps").at("evcount")->Fill("req. 1 proton",1);



    //fitToEtaPrime.SetEgammaBeam(data.TaggerHits().at(0)->PhotonEnergy());
//    fitToEtaPrime.SetProtonTAPS(protons.at(0));
    //fitToEtaPrime.SetProtonTAPS(mcproton);
    //fitToEtaPrime.SetPhotons(photons);

    //result_fitToEtaPrime = fitToEtaPrime.DoFit();


    result_t result_3pi0    = Make3pi0(photons);
    result_t result_eta2pi0 = MakeEta2pi0(photons);
    result_t result_mc      = MakeMC3pi0(mcdata);



    const double chi2cut(3);

    if (result_3pi0.chi2() < chi2cut )
    {
        FillIm(result_3pi0, ParticleTypeDatabase::Pi0, (TH1D*) hists.at("signal").at("IM_pi0"));
        FillImEtaPrime(result_3pi0,(TH1D*)hists.at("signal").at("IM_etap"));
        hists.at("channels").at("signal_chi2")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree()).c_str(),1);
        hists.at("steps").at("evcount")->Fill("#chi^{2} cut signal",1);

        /*
        hists.at("kinfit").at("signal_chi2")->Fill(result_fitToEtaPrime.ChiSquare);
        hists.at("kinfit").at("signal_niter")->Fill(result_fitToEtaPrime.NIterations);
        auto& kinfitvars = result_fitToEtaPrime.Variables;
        hists.at("kinfit").at("signal_egamma_before")->Fill(kinfitvars.at(fitToEtaPrime.egammaName).Value.Before);
        hists.at("kinfit").at("signal_egamma_after")->Fill(kinfitvars.at(fitToEtaPrime.egammaName).Value.After);
        */
        for (const auto& ph: photons)
            hists["P"]["gamma_signal"]->Fill(ph->Theta()*TMath::RadToDeg(),ph->Ek());
    }

    if (result_eta2pi0.chi2() < chi2cut)
    {
        FillIm(result_eta2pi0, ParticleTypeDatabase::Pi0,(TH1D*) hists.at("ref").at("IM_pions"));
        FillIm(result_eta2pi0, ParticleTypeDatabase::Eta,(TH1D*) hists.at("ref").at("IM_etas"));
        FillImEtaPrime(result_eta2pi0,(TH1D*) hists.at("ref").at("IM_etap"));
        hists.at("channels").at("ref_chi2")->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree()).c_str(),1);
        hists.at("steps").at("evcount")->Fill("#chi^{2} cut reference",1);
        for (const auto& ph: photons)
            hists["P"]["gamma_ref"]->Fill(ph->Theta()*TMath::RadToDeg(),ph->Ek());
    }

    if ( result_mc.success)
    {
        DalitzVars channel(result_mc);
        hists.at("dalitz").at("mc_xy")->Fill(channel.s1,channel.s3);
        hists.at("dalitz").at("mc_z")->Fill(channel.z);
    }
    if (result_3pi0.chi2() < chi2cut)
    {
        DalitzVars channel(result_3pi0);
        hists.at("dalitz").at("xy")->Fill(channel.s1,channel.s3);
        hists.at("dalitz").at("z")->Fill(channel.z);
    }

}


void Etap3pi0::ShowResult()
{
    for (auto& category: hists)
    {
        canvas c(category.first);
        for (auto& h: category.second)
            c << h.second;
        c << endc;
    }

    canvas("Dalitz-Plots")      << drawoption("colz") << hists.at("dalitz").at("xy") << hists.at("dalitz").at("mc_xy") << endc;
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

    TMean = 0;
    for (const auto& meson: r.mesons)
        TMean += meson.first->Ek();
    TMean = TMean / 3.0;

    x = ( r.mesons[0].first->Ek() - r.mesons[1].first->Ek() ) / ( TMath::Sqrt(3) * TMean);

//    y = 0;
//    for (const auto& meson: r.mesons)
//        y += meson.first->M();
//    y = y / ( 3 * r.mother.M());
    y = 1;
    y = (y * r.mesons[2].first->Ek() / TMean) - 1;

    z = std_ext::sqr(x) + std_ext::sqr(y);

}



double Etap3pi0::KinFitter::kinVector::energySmear(const double& E) const
{
    return  0.02 * E * std::pow(E,-0.36);
}

void Etap3pi0::KinFitter::kinVector::SetEkThetaPhi(double ek, double theta, double phi)
{
    Ek     = ek;
    Theta  = theta;
    Phi    = phi;
    sEk    = energySmear(ek);
    sTheta = 2.5 * TMath::DegToRad();
    if ( Theta > 20 * TMath::DegToRad() && Theta < 160 * TMath::DegToRad())
    {
        sPhi = sTheta / std::sin(Theta);
    }
    else
    {
        sPhi = 1 * TMath::DegToRad();
    }
}

TLorentzVector Etap3pi0::KinFitter::GetVector(const vector<double>& EkThetaPhi, const double m) const
{
    TLorentzVector lv;
    const double E = sqrt(std_ext::sqr(EkThetaPhi[0] + m) - std_ext::sqr(m));

    lv.SetE(E);
    lv.SetTheta(EkThetaPhi[1]);
    lv.SetPhi(EkThetaPhi[2]);

    return lv;
}

void Etap3pi0::KinFitter::SetEgammaBeam(const double& ebeam)
{
    EgammaBeam.first = ebeam;
    EgammaBeam.second = taggerSmear(ebeam);
}

void Etap3pi0::KinFitter::SetProtonTAPS(const data::ParticlePtr& proton)
{
    ProtonTAPS.SetEkThetaPhi(proton->Ek(),proton->Theta(),proton->Phi());
}

void Etap3pi0::KinFitter::SetPhotons(const std::vector<ParticlePtr>& photons)
{
    assert(Photons.size() == photons.size());

    for ( unsigned i = 0 ; i < Photons.size() ; ++ i)
        Photons.at(i).SetEkThetaPhi(
                    photons.at(i)->Ek(),
                    photons.at(i)->Theta(),
                    photons.at(i)->Phi());
}

double Etap3pi0::KinFitter::taggerSmear(const double& E) const
{
    return  0.02 * E * std::pow(E,-0.36);
}

void Etap3pi0::AddHist1D(
        const std::string& category, const std::string& hname,
        const std::string& title,
        const std::string& xlabel, const std::string& ylabel,
        const BinSettings& bins
        )
{
    hists[category][hname] = HistFac.makeTH1D(title,xlabel,ylabel,bins,(category + string("_") + hname));
}

void Etap3pi0::AddHist2D(const string& category, const string& hname,
                         const string& title,
                         const string& xlabel, const string& ylabel,
                         const BinSettings& xbins, const BinSettings& ybins)
{
    hists[category][hname] = HistFac.makeTH2D(title,xlabel,ylabel,xbins,ybins,(category + string("_") + hname));
}

Etap3pi0::KinFitter::KinFitter(const ParticleTypeDatabase::Type& motherParticle):
    aplcon("Kinfitter_photons"),
    ProtonTAPS("ProtonTAPS"),
    Photons({kinVector("Photon0"),
             kinVector("Photon1"),
             kinVector("Photon2"),
             kinVector("Photon3"),
             kinVector("Photon4"),
             kinVector("Photon5")}),
    IM_Mother(motherParticle.Mass())
{

    aplcon.LinkVariable(egammaName,
                        { addressof(EgammaBeam.first) },
                        { addressof(EgammaBeam.second)});
    aplcon.LinkVariable(ProtonTAPS.Name,
                        ProtonTAPS.Adresses(),
                        ProtonTAPS.Adresses_Sigma());

    vector<string> namesLInv      = { egammaName, ProtonTAPS.Name };
    vector<string> namesEtapMass;

    for ( auto& photon: Photons)
    {
        aplcon.LinkVariable(photon.Name,
                            photon.Adresses(),
                            photon.Adresses_Sigma());
        namesLInv.push_back(photon.Name);
        namesEtapMass.push_back(photon.Name);
    }

    auto LorentzInvariance = [this] (const vector<vector<double>> values)
    {
        //  Beam-LV:
        TLorentzVector constraint(0,0,values.at(0)[0],ParticleTypeDatabase::Proton.Mass());
        constraint = constraint - this->GetVector(values[1],ParticleTypeDatabase::Proton.Mass());
        for ( unsigned i = 0 ; i < 6 ; ++ i)
            constraint = constraint - this->GetVector(values[i + 2],0);
        return vector<double>(
               { constraint.X(),
                 constraint.Y(),
                 constraint.Z(),
                 constraint.Z()  } );
    };

    aplcon.AddConstraint("LInv",namesLInv,LorentzInvariance);

    auto etapMass = [this] ( const vector<vector<double>> values )
    {
        TLorentzVector SumPhotons(0,0,0,0);
        for ( unsigned i = 0 ; i < 6 ; ++ i)
            SumPhotons = SumPhotons + this->GetVector(values[i],0);

        return SumPhotons.M2() - std_ext::sqr(IM_Mother);

    };

    aplcon.AddConstraint("etapMass",namesEtapMass,etapMass);

}

AUTO_REGISTER_PHYSICS(Etap3pi0)
