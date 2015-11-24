#include "kinFit.h"
#include <APLCON.hpp>
#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "TVector3.h"
#include <cassert>
#include "base/ParticleType.h"

using namespace std;
using namespace ant::std_ext;
using namespace ant;



double KinFitter::EnergyResolution(const analysis::data::ParticlePtr& p) const
{
    if(p->Candidate()) {
        if(p->Candidate()->Detector() & Detector_t::Type_t::CB) {

            return fct_GlobalEResolution(p->Ek());

        } if(p->Candidate()->Detector() & Detector_t::Type_t::TAPS) {

            return fct_GlobalEResolution(p->Ek()); ///@todo check TAPS energy resolution fct
        }
    } else {
        throw runtime_error("Particle without candidate!");
    }
    return 0.0;
}

double KinFitter::ThetaResolution(const analysis::data::ParticlePtr&) const
{
    return degree_to_radian(2.5);
}

double KinFitter::PhiResolution(const analysis::data::ParticlePtr& p) const
{
    if(p->Candidate()) {
        if(p->Candidate()->Detector() & Detector_t::Type_t::CB) {

            return p->Candidate()->Theta() / std::sin(p->Candidate()->Theta());

        } if(p->Candidate()->Detector() & Detector_t::Type_t::TAPS) {   ///@todo check TAPS Theta resolution

            return p->Candidate()->Theta() / std::sin(p->Candidate()->Theta());
        }
    } else {
        throw runtime_error("Particle without candidate!");
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
    return  0.02 * E * std::pow(E,-0.36);
}

double KinFitter::fct_TaggerEGausSigma(double)
{
    return  3.0/sqrt(12.0);
}

KinFitter::KinFitter():
    aplcon(make_unique<APLCON>("OmegaFitter"))
{

    aplcon->LinkVariable(Beam.name,
                        { Beam.Adresses() },
                        { Beam.Adresses_Sigma()});

    aplcon->LinkVariable(Proton.Name,
                        Proton.Adresses(),
                        Proton.Adresses_Sigma());

    vector<string> namesLInv      = { Beam.name, Proton.Name };

    for ( auto& photon: Photons)
    {
        aplcon->LinkVariable(photon.Name,
                            photon.Adresses(),
                            photon.Adresses_Sigma());

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

        for ( unsigned photon = 0 ; photon < 3 ; ++ photon)
            constraint -= KinFitter::GetVector(values[photon + 2], 0.0);

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
    Proton.sEk   = 0.0; // unmeasured

    Proton.Theta  = proton->Theta();
//    Proton.sTheta = ...

    Proton.Phi   = proton->Phi();
    Proton.sPhi  = degree_to_radian(1.0);

}

void KinFitter::SetPhotons(const std::vector<analysis::data::ParticlePtr> &photons_data)
{
    assert(Photons.size() == photons_data.size());

    for ( unsigned i = 0 ; i < Photons.size() ; ++ i) {
        auto& photon = Photons.at(i);
        auto& data   = photons_data.at(i);

        photon.Ek  = data->Ek();
        photon.sEk = EnergyResolution(data);

        photon.Theta  = data->Theta();
        photon.sTheta = ThetaResolution(data);

        photon.Phi  = data->Phi();
        photon.sPhi = PhiResolution(data);
    }
}

APLCON::Result_t KinFitter::DoFit() { return aplcon->DoFit(); }


KinFitter::PhotonBeamVector::PhotonBeamVector(const string& Name):
    name(Name)
{

}
