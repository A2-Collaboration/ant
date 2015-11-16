#include "kinFit.h"
#include <APLCON.hpp>
#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "TVector3.h"
#include <cassert>

using namespace std;
using namespace ant::std_ext;
using namespace ant;



TLorentzVector KinFitter::GetVector(const std::vector<double>& EkThetaPhi, const double m)
{
    const mev_t E = EkThetaPhi[0] + m;
    const mev_t p = sqrt( sqr(E) - sqr(m) );

    /// \bug This might be inefficient...

    TVector3 pv(1,0,0);

    pv.SetMagThetaPhi(p,EkThetaPhi[1],EkThetaPhi[2]);

    return TLorentzVector(pv,E);
}

KinFitter::KinFitter():
    aplcon(make_unique<APLCON>("OmegaFitter"))
{

    aplcon->LinkVariable(egammaName,
                        { addressof(EgammaBeam.first) },
                        { addressof(EgammaBeam.second)});

    aplcon->LinkVariable(Proton.Name,
                        Proton.Adresses(),
                        Proton.Adresses_Sigma());

    vector<string> namesLInv      = { egammaName, Proton.Name };

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
        TLorentzVector constraint(0, 0, Ebeam, ParticleTypeDatabase::Proton.Mass() + Ebeam);

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
    EgammaBeam.first = ebeam;
    EgammaBeam.second = taggerSmear(ebeam);
}

void KinFitter::SetProton(const analysis::data::ParticlePtr &proton)
{
    Proton.SetEkThetaPhi(proton->Ek(),proton->Theta(),proton->Phi());
}

void KinFitter::SetPhotons(const std::vector<analysis::data::ParticlePtr> &photons)
{
    assert(Photons.size() == photons.size());

    for ( unsigned i = 0 ; i < Photons.size() ; ++ i)
        Photons.at(i).SetEkThetaPhi(
                    photons.at(i)->Ek(),
                    photons.at(i)->Theta(),
                    photons.at(i)->Phi());
}

APLCON::Result_t KinFitter::DoFit() { return aplcon->DoFit(); }

double KinFitter::kinVector::energySmear(const double &E) const
{
    return  0.02 * E * std::pow(E,-0.36);
}

double KinFitter::taggerSmear(const double &E) const
{
    return  0.02 * E * std::pow(E,-0.36);
}

void KinFitter::kinVector::SetEkThetaPhi(double ek, double theta, double phi)
{
    this->Ek = ek;
    this->Theta = theta;
    this->Phi = phi;

    this->sEk = energySmear(ek);
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
