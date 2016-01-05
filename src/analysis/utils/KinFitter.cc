#include "KinFitter.h"

#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "base/ParticleType.h"

#include "base/WrapTFile.h"
#include <cassert>
#include "expconfig/ExpConfig.h"

#include <APLCON.hpp>

#include "TTree.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TF1.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;


double KinFitter::EnergyResolution(const analysis::data::ParticlePtr& p) const
{
    if(p->Candidate) {
        if(p->Candidate->Detector & Detector_t::Type_t::CB) {

            return 1.07134e-02 * p->Ek();

        } if(p->Candidate->Detector & Detector_t::Type_t::TAPS) {

            return 3.5E-2 * p->Ek();
        }
    } else {
        if(p->Theta() < degree_to_radian(20.0) ) {

            return 1.07134e-02 * p->Ek();

        } else {

            return 3.5E-2 * p->Ek();
        }
    }
    return 0.0;
}

double KinFitter::ThetaResolution(const analysis::data::ParticlePtr& p) const
{
    if(p->Candidate) {
        if(p->Candidate->Detector & Detector_t::Type_t::CB) {

            return cb_sigma_theta.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());

        } if(p->Candidate->Detector & Detector_t::Type_t::TAPS) {   ///@todo check TAPS Theta resolution
             return degree_to_radian(2.5);
            //return taps_sigma_theta.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());
        }
    } else {
        return degree_to_radian(2.5);
    }
    return 0.0;
}

double KinFitter::PhiResolution(const analysis::data::ParticlePtr& p) const
{
    if(p->Candidate) {
        if(p->Candidate->Detector & Detector_t::Type_t::CB) {

            return cb_sigma_phi.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());

        } if(p->Candidate->Detector & Detector_t::Type_t::TAPS) {   ///@todo check TAPS Theta resolution
            return p->Theta() / std::sin(p->Theta());
            //return taps_sigma_phi.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());
        }
    } else {
        if(p->Theta() < degree_to_radian(20.0) ) {

            return p->Theta() / std::sin(p->Theta());

        } else {   ///@todo check TAPS Theta resolution

            return p->Theta() / std::sin(p->Theta());
        }
    }
    return 0.0;
}

TLorentzVector KinFitter::GetVector(const std::vector<double>& EkThetaPhi, const double m)
{
    const mev_t E = EkThetaPhi[0] + m;
    const mev_t p = sqrt( sqr(E) - sqr(m) );

    /// \bug This might be inefficient...

    TVector3 pv(1,0,0);

    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);

    return TLorentzVector(pv,E);
}

double KinFitter::fct_GlobalEResolution(double E)
{
 //   return  0.02 * E * std::pow(E,-0.36);
    return 1.07134e-02 * E; // roughly extracted
}

double KinFitter::fct_TaggerEGausSigma(double)
{
    return  3.0/sqrt(12.0);
}

KinFitter::KinFitter(const std::string& name, unsigned numGammas):
    aplcon(make_unique<APLCON>(name))
{

    Photons.reserve(numGammas);
    for(unsigned i=0; i<numGammas;++i) {
        Photons.emplace_back(kinVector("Photon"+to_string(i)));
    }

    aplcon->LinkVariable(Beam.name,
                        { Beam.Adresses() },
                        { Beam.Adresses_Sigma() });

    aplcon->LinkVariable(Proton.Name,
                         Proton.Addresses(),
                         Proton.Addresses_Sigma(),
                         Proton.Addresses_Pulls());

    vector<string> namesLInv      = { Beam.name, Proton.Name };

    for ( auto& photon: Photons)
    {
        aplcon->LinkVariable(photon.Name,
                             photon.Addresses(),
                             photon.Addresses_Sigma(),
                             photon.Addresses_Pulls());

        namesLInv.push_back(photon.Name);
    }

    auto LorentzInvariance = [] (const vector<vector<double>> values)
    {
        // beam    TLorentzVector(0.0, 0.0, PhotonEnergy(), PhotonEnergy());
        // target  TLorentzVector(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
        const auto& Ebeam  = values[0][0];
        const auto& proton = values[1];

        //  Beam-LV:
        const TLorentzVector beam(0, 0, Ebeam, Ebeam);
        const TLorentzVector tg(0,0,0,ParticleTypeDatabase::Proton.Mass());

        TLorentzVector constraint = tg + beam;

        constraint -= KinFitter::GetVector(proton, ParticleTypeDatabase::Proton.Mass());

        const auto s = values.size();
        for ( unsigned photon = 0 ; photon < s-2 ; ++ photon)
            constraint -= KinFitter::GetVector(values[photon + 2], ParticleTypeDatabase::Photon.Mass());

        return vector<double>(
               { constraint.X(),
                 constraint.Y(),
                 constraint.Z(),
                 constraint.E()} );
    };

    aplcon->AddConstraint("LInv",{namesLInv},LorentzInvariance);

}

void KinFitter::SetEgammaBeam(const double &ebeam)
{
    Beam.E     = ebeam;
    Beam.Sigma = fct_TaggerEGausSigma(ebeam);
}

void KinFitter::SetProton(const analysis::data::ParticlePtr &proton)
{
    Proton.Ek    = proton->Ek();
    Proton.sigmaEk   = 0.0; // unmeasured

    Proton.Theta  = proton->Theta();
    Proton.sigmaTheta = ThetaResolution(proton)*5.0;

    Proton.Phi   = proton->Phi();
    Proton.sigmaPhi  = PhiResolution(proton)*5.0;

}

data::ParticlePtr KinFitter::GetFittedProton() const
{
    return make_shared<data::Particle>(ParticleTypeDatabase::Proton, Proton.Ek, Proton.Theta, Proton.Phi);
}

data::ParticleList KinFitter::GetFittedPhotons() const
{
    data::ParticleList photons;
    for(const auto& photon : Photons) {
        photons.emplace_back(make_shared<data::Particle>(ParticleTypeDatabase::Photon,
                                                         photon.Ek, photon.Theta, photon.Phi));
    }
    return photons;
}



void KinFitter::SetPhotons(const std::vector<analysis::data::ParticlePtr> &photons_data)
{
    assert(Photons.size() == photons_data.size());

    for ( unsigned i = 0 ; i < Photons.size() ; ++ i) {
        auto& photon = Photons.at(i);
        auto& data   = photons_data.at(i);

        photon.Ek  = data->Ek();
        photon.sigmaEk = EnergyResolution(data);

        photon.Theta  = data->Theta();
        photon.sigmaTheta = ThetaResolution(data);

        photon.Phi  = data->Phi();
        photon.sigmaPhi = PhiResolution(data);
    }
}

APLCON::Result_t KinFitter::DoFit() {

    const auto res = aplcon->DoFit();

    result_chi2ndof    = res.ChiSquare / res.NDoF;
    result_iterations  = res.NIterations;
    result_status      = static_cast<int>(res.Status);
    result_probability = res.Probability;

    return res;
}

void KinFitter::SetupBranches(TTree* tree, string branch_prefix)
{
    if(branch_prefix.empty())
        branch_prefix = aplcon->GetName();

    Proton.SetupBranches(tree, branch_prefix);
    for(auto& p : Photons) {
        p.SetupBranches(tree, branch_prefix);
    }

    tree->Branch((branch_prefix+"_chi2dof").c_str(),     &result_chi2ndof);
    tree->Branch((branch_prefix+"_iterations").c_str(),  &result_iterations);
    tree->Branch((branch_prefix+"_status").c_str(),      &result_status);
    tree->Branch((branch_prefix+"_probability").c_str(), &result_probability);
}

void KinFitter::LoadSigmaData(const string& filename)
{
    const auto& setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw std::runtime_error("No Setup found!");


    WrapTFileInput f(filename);

    const auto& CB = setup->GetDetector(Detector_t::Type_t::CB);
    if(!CB)
        throw std::runtime_error("No CB detector defined in setup");

    const auto nCB = CB->GetNChannels();

    const auto& TAPS = setup->GetDetector(Detector_t::Type_t::TAPS);
    if(!TAPS)
        throw std::runtime_error("No TAPS detector defined in setup");

    const auto nTAPS = TAPS->GetNChannels();

    cb_sigma_theta.Load(f, "cb_sigma_theta", int(nCB));
    cb_sigma_phi.Load(  f, "cb_sigma_phi",   int(nCB));

//    taps_sigma_theta.Load(f, "taps_sigma_theta", int(nTAPS));
//    taps_sigma_phi.Load(  f, "taps_sigma_phi",   int(nTAPS));
}


KinFitter::PhotonBeamVector::PhotonBeamVector(const string& Name):
    name(Name)
{

}

void KinFitter::kinVector::SetupBranches(TTree* tree, const string& prefix)
{
    tree->Branch((prefix+"_"+Name+"_Ek_pull").c_str(),     addressof(pullEk));
    tree->Branch((prefix+"_"+Name+"_Theta_pull").c_str(),  addressof(pullTheta));
    tree->Branch((prefix+"_"+Name+"_Phi_pull").c_str(),    addressof(pullPhi));
    tree->Branch((prefix+"_"+Name+"_Ek_sigma").c_str(),    addressof(sigmaEk));
    tree->Branch((prefix+"_"+Name+"_Theta_sigma").c_str(), addressof(sigmaTheta));
    tree->Branch((prefix+"_"+Name+"_Phi_sigma").c_str(),   addressof(sigmaPhi));

}

double KinFitter::angular_sigma::GetSigma(const unsigned element, const double E) const
{
    const double vp0 = p0->GetBinContent(int(element));
    const double vp1 = p1->GetBinContent(int(element));
    const double vp2 = p2->GetBinContent(int(element));

    return f(E, vp0, vp1, vp2);

}

double KinFitter::angular_sigma::f(const double x, const double p0, const double p1, const double p2) noexcept
{
    //return p0*exp(p1*x)+p2;
    return exp(p0 + p1*x) + p2;
}

double KinFitter::angular_sigma::f_root(const double* x, const double* p) noexcept
{
    return f(x[0], p[0], p[1], p[2]);
}

TF1* KinFitter::angular_sigma::GetTF1(const std::string& name)
{
    return new TF1(name.c_str(), f_root, 0, 1600, 3);
}

KinFitter::angular_sigma::~angular_sigma()
{}

void KinFitter::angular_sigma::Load(WrapTFile& f, const string& prefix, const int bins)
{
    p0 = LoadHist(f, prefix+"_p0", bins);
    p1 = LoadHist(f, prefix+"_p1", bins);
    p2 = LoadHist(f, prefix+"_p2", bins);
}

KinFitter::angular_sigma::Hist KinFitter::angular_sigma::LoadHist(WrapTFile& f, const string& name, const int bins)
{
    auto h = f.GetSharedHist<TH1D>(name);

    if(!h)
        throw std::runtime_error("Could not load histogram " + name);

    if(h->GetNbinsX() != bins)
        throw std::runtime_error(formatter() << "TH1D " << name << " has wrong number of bins: " << h->GetNbinsX()  << ", expected: " << bins);

    return h;
}

KinFitter::angular_sigma::angular_sigma()
{

}
