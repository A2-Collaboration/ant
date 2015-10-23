#include "TestAPLCON.h"
#include "data/Particle.h"
#include "data/Candidate.h"
#include "plot/root_draw.h"
#include <string>
#include "utils/combinatorics.h"
#include <vector>
#include <numeric>
#include <functional>
#include <APLCON.hpp>
#include <iomanip>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

std::default_random_engine TestAPLCON::FitParticle::generator;

TLorentzVector TestAPLCON::FitParticle::Make(const std::vector<double> &EkThetaPhi, const Double_t m) {
    const double E = EkThetaPhi[0] + m;
    const Double_t p = sqrt( E*E - m*m );
    TVector3 pv(1,0,0);
    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
    TLorentzVector l(pv, E);
    return l;
}

void TestAPLCON::FitParticle::Smear() {
    // set the sigmas here,
    // then the fitter knows them as well (because they're linked)

    Ek_Sigma = 0.02*Ek*pow(Ek,-0.36);
    Theta_Sigma = 2.5*TMath::DegToRad();
    if(Theta>20*TMath::DegToRad() && Theta<160*TMath::DegToRad()) {
        Phi_Sigma = Theta_Sigma/sin(Theta);
    }
    else {
        Phi_Sigma = 1*TMath::DegToRad();
    }

    // then artificially smear the values with gaussians
    using gauss_t = std::normal_distribution<double>;
    gauss_t gauss_Ek(0, Ek_Sigma);
    Ek += gauss_Ek(generator);
    gauss_t gauss_Theta(0, Theta_Sigma);
    Theta += gauss_Theta(generator);
    gauss_t gauss_Phi(0, Phi_Sigma);
    Phi += gauss_Phi(generator);
}

void TestAPLCON::FillIM(TH1D *h, const std::vector<TestAPLCON::FitParticle> &photons) {
    TLorentzVector sum(0,0,0,0);
    for(const auto& p : photons) {
        sum += FitParticle::Make(p, ParticleTypeDatabase::Photon.Mass());
    }
    h->Fill(sum.M());
}

TestAPLCON::TestAPLCON(const string& name, PhysOptPtr opts) :
    Physics(name, opts),
    fitter(name),
    photons(nPhotons)
{


    const BinSettings energy_bins(1000,0,1600);
    const BinSettings tagger_bins(2000,0.0,2000);
    const BinSettings ntaggerhits_bins(100);
    const BinSettings veto_bins(1000,0,10.0);
    const BinSettings particle_bins(10,0,10);
    const BinSettings particlecount_bins(16,0,16);
    const BinSettings chisqare_bins(100,0,30);
    const BinSettings probability_bins(100,0,1);
    const BinSettings iterations_bins(15,0,15);
    const BinSettings im_bins(200,IM-100,IM+100);
    const BinSettings vertex_bins(200,-10,10);

    banana = HistFac.makeTH2D(
                "PID Bananas",
                "Cluster Energy [MeV]",
                "Veto Energy [MeV]",
                energy_bins,
                veto_bins,
                "pid"
                );

    particles = HistFac.makeTH1D(
                "Identified particles",
                "Particle Type",
                "#",
                particle_bins,
                "ParticleTypes"
                );
    tagger = HistFac.makeTH1D(
                "Tagger Spectrum",
                "Photon Beam Energy",
                "#",
                tagger_bins,
                "TaggerSpectrum"
                );

    ntagged = HistFac.makeTH1D(
                "Tagger Hits",
                "Tagger Hits / event",
                "#",
                ntaggerhits_bins,
                "nTagged"
                );

    cbesum = HistFac.makeTH1D(
                "CB Energy Sum",
                "E [MeV]",
                "#",
                energy_bins,
                "esum"
                );

    for( auto& t : ParticleTypeDatabase::DetectableTypes() ) {
        numParticleType[t]= HistFac.makeTH1D("Number of "+t->PrintName(),
                                      "number of "+t->PrintName()+"/ event",
                                      "",
                                      particlecount_bins);
    }

    // setup fitter for nPhotons

    fitter.LinkVariable("Beam",    beam.Link(),       beam.LinkSigma());
    fitter.LinkVariable("Proton",  proton.Link(),     proton.LinkSigma());

    vector<string> photon_names;
    for(size_t i=0;i<nPhotons;i++) {
        stringstream s_photon;
        s_photon << "Photon" << (i+1);
        photon_names.push_back(s_photon.str());
        fitter.LinkVariable(s_photon.str(), photons[i].Link(), photons[i].LinkSigma());
    }
    vector<string> all_names = {"Beam", "Proton"};
    all_names.insert(all_names.end(),photon_names.begin(),photon_names.end());

    // Constraint: Incoming 4-vector = Outgoing 4-vector
    auto EnergyMomentumBalance = [] (const vector< vector<double> >& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, ParticleTypeDatabase::Proton.Mass());
        // assume first particle is beam photon
        TLorentzVector diff = target + FitParticle::Make(particles[0], ParticleTypeDatabase::Photon.Mass());
        // assume second particle outgoing proton
        diff -= FitParticle::Make(particles[1], ParticleTypeDatabase::Proton.Mass());
        // subtract the rest, assumed to be photons
        for(size_t i=2;i<particles.size();i++) {
            diff -= FitParticle::Make(particles[i], ParticleTypeDatabase::Photon.Mass());
        }

        return {diff.X(), diff.Y(), diff.Z(), diff.T()};

    };
    fitter.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification
    auto RequireIM = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(const auto& p : photons) {
            sum += FitParticle::Make(p, ParticleTypeDatabase::Photon.Mass());
        }
        return sum.M() - IM;
    };
    if(includeIMconstraint)
        fitter.AddConstraint("RequireIM", photon_names, RequireIM);

    // Constraint: Vertex position in z direction: v_z (positive if upstream)
    // if the photon originated from (0,0,v_z) instead of origin,
    // the corrected angle theta' is given by
    // tan(theta') = (R sin(theta))/(R cos(theta) - v_z)
    // R is the CB radius, 10in aka 25.4cm

    auto VertexConstraint = [&] (vector< vector<double> >& photons) -> double
    {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = photons.back()[0];
        photons.resize(photons.size()-1); // get rid of last element
        // correct each photon's theta angle,
        // then calculate invariant mass of all photons
        TLorentzVector sum(0,0,0,0);
        for(auto& p : photons) {
            const double theta = p[1]; // second element is theta
            const double theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
            p[1] = theta_p;
            sum += FitParticle::Make(p, ParticleTypeDatabase::Photon.Mass());
        }
        return sum.M() - IM;
    };

    if(includeVertexFit) {
        fitter.AddUnmeasuredVariable("v_z"); // default value 0
        fitter.AddConstraint("VertexConstraint", photon_names + std::vector<string>{"v_z"}, VertexConstraint);
    }

    static_assert(!(includeIMconstraint && includeVertexFit), "Do not enable Vertex and IM Fit at the same time");

    // make fitter histograms
    chisquare   = HistFac.makeTH1D("ChiSqare","ChiSquare","#",chisqare_bins,"chisquare");
    probability = HistFac.makeTH1D("Probability","Probability","#",probability_bins,"probability");
    iterations = HistFac.makeTH1D("Number of iterations","Iterations","#",iterations_bins,"iterations");

    // pull histograms are created on the fly

    stringstream ng;
    ng << nPhotons << "g";
    im_true = HistFac.makeTH1D("IM "+ng.str()+" true","IM","#",im_bins,"im_true");
    im_smeared = HistFac.makeTH1D("IM "+ng.str()+" smeared","IM","#",im_bins,"im_smeared");
    im_fit = HistFac.makeTH1D("IM "+ng.str()+" fit","IM","#",im_bins,"im_fit");

    vertex_z_before =  HistFac.makeTH1D("Vertex Z Before","v_z / cm","#",vertex_bins,"vertex_z_before");
    vertex_z_after =  HistFac.makeTH1D("Vertex Z After","v_z / cm","#",vertex_bins,"vertex_z_after");

    APLCON::Fit_Settings_t settings = fitter.GetSettings();
    settings.MaxIterations = 50;
    fitter.SetSettings(settings);

    cout.precision(3);
    APLCON::PrintFormatting::Width = 11;
}


void TestAPLCON::ProcessEvent(const Event &event)
{


    for(auto& cand : event.Reconstructed().Candidates()) {
        banana->Fill(cand->ClusterEnergy(), cand->VetoEnergy());
    }

    for(auto& particle : event.Reconstructed().Particles().GetAll()) {
        particles->Fill(particle->Type().PrintName().c_str(), 1);
    }

    ntagged->Fill(event.Reconstructed().TaggerHits().size());

    cbesum->Fill(event.Reconstructed().TriggerInfos().CBEenergySum());

    for( auto& t : ParticleTypeDatabase::DetectableTypes() ) {
        try {
            numParticleType.at(t)->Fill(event.Reconstructed().Particles().Get(*t).size());
        } catch (...) {}
    }

    for(const auto& taggerhit : event.MCTrue().TaggerHits()) {
        tagger->Fill(taggerhit->PhotonEnergy());

        // find the photons and one proton
        size_t foundPhotons = 0;
        for(const auto& p : event.MCTrue().Particles().GetAll()) {
            if(p->Type() == ParticleTypeDatabase::Proton) {
                proton.SetFromVector(*p);
            }
            else if(foundPhotons<nPhotons && p->Type() == ParticleTypeDatabase::Photon) {
                photons[foundPhotons].SetFromVector(*p);
                foundPhotons++;
           }
        }
        if(foundPhotons != nPhotons)
            continue;

        beam.SetFromVector(taggerhit->PhotonBeam());

        FillIM(im_true, photons);

        // smear the MC true data
        proton.Smear();
        for(auto& photon : photons)
            photon.Smear();
        beam.Smear();

        FillIM(im_smeared, photons);

        // let APLCON do the work
        const APLCON::Result_t& result = fitter.DoFit();

        //cout << result << endl;

        if(result.Status != APLCON::Result_Status_t::Success) {
            //cout << result << endl;
            continue;
        }

        for(const auto& it_map : result.Variables) {
            const string& varname = it_map.first;
            const APLCON::Result_Variable_t& var = it_map.second;
            auto it_pull = pulls.find(varname);
            TH1D* h_pull;
            if(it_pull == pulls.end()) {
                // not found so far, create the histogram on-the-fly
                const BinSettings pull_bins(50,-3,3);
                stringstream title;
                title << "Pull " << var.PristineName << " " << component.at(var.Index);
                h_pull = HistFac.makeTH1D(title.str(),
                                           "Pull", "#",
                                           pull_bins,
                                           "pull_"+varname);
                pulls[varname] = h_pull;
            }
            else {
                h_pull = it_pull->second;
            }
            h_pull->Fill(var.Pull);
        }
        chisquare->Fill(result.ChiSquare);
        probability->Fill(result.Probability);
        iterations->Fill(result.NIterations);

        if(includeVertexFit) {
            vertex_z_after->Fill(result.Variables.at("v_z").Value.After);
            vertex_z_before->Fill(result.Variables.at("v_z").Value.Before);
        }

        FillIM(im_fit, photons);
    }

}

void TestAPLCON::Finish()
{

}

void TestAPLCON::ShowResult()
{
    canvas c("TestAPLCON: Overview");
    c << drawoption("colz") << banana
      << padoption::set(padoption_t::Legend)
      << particles
      << padoption::unset(padoption_t::Legend)
      << tagger << ntagged << cbesum << endc;

    canvas c_pulls("TestAPLCON: Pulls");
    c_pulls << padoption::set(padoption_t::LogY);
    for( auto& p : pulls ) {
        c_pulls << p.second;
    }
    c_pulls << endc;


    canvas c_fitter("TestAPLCON: Fitter");
    c_fitter << chisquare << probability << iterations
             << im_true << im_smeared << im_fit
             << vertex_z_before << vertex_z_after << endc;

}



AUTO_REGISTER_PHYSICS(TestAPLCON)
