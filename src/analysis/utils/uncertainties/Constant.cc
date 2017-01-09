#include "Constant.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;
using namespace ant::analysis::utils::UncertaintyModels;

Constant::Constant()
{}

Constant::~Constant()
{}

std::shared_ptr<Constant> Constant::make()
{
    return std::make_shared<Constant>();
}

Uncertainties_t Constant::GetSigmas(const TParticle &particle) const
{
    if(particle.Candidate->Detector & Detector_t::Type_t::CB) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            return photon_cb;
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            return proton_cb;
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

    } else if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            return photon_taps;
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            return proton_taps;
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }
    }
    else {
        throw Exception("Unexpected Detector: " + string(particle.Candidate->Detector));
    }
}

ConstantRelativeE::ConstantRelativeE()
{}

ConstantRelativeE::~ConstantRelativeE()
{}

Uncertainties_t ConstantRelativeE::GetSigmas(const TParticle &particle) const
{
    auto s = Constant::GetSigmas(particle);
    s.sigmaEk *= particle.Ek();

    return s;
}

std::shared_ptr<ConstantRelativeE> ConstantRelativeE::makeMCLongTarget()
{
    auto s = std::make_shared<ConstantRelativeE>();

    s->photon_cb   = { 0.0107, std_ext::degree_to_radian(3.79), std_ext::degree_to_radian(1.78)};
    s->photon_taps = { 0.035,  std_ext::degree_to_radian(0.42), std_ext::degree_to_radian(1.15)};

    s->proton_cb   = { 0.0,    std_ext::degree_to_radian(5.5), std_ext::degree_to_radian(5.3)};
    s->proton_taps = { 0.0,    std_ext::degree_to_radian(2.8), std_ext::degree_to_radian(4.45)};

    return s;
}

ConstantRelativeEpow::ConstantRelativeEpow()
{}

ConstantRelativeEpow::~ConstantRelativeEpow()
{}

Uncertainties_t ConstantRelativeEpow::GetSigmas(const TParticle &particle) const
{
    Uncertainties_t s;

    if(particle.Candidate->Detector & Detector_t::Type_t::CB) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            s = photon_cb;
            s.sigmaEk = s.sigmaEk* particle.Ek() * pow(particle.Ek(), Eexp_cb);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            s = proton_cb;
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

    } else if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            s = photon_taps;
            s.sigmaEk = s.sigmaEk* particle.Ek() * pow(particle.Ek(), Eexp_taps);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            s = proton_taps;
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }
    }
    else {
        throw Exception("Unexpected Detector: " + string(particle.Candidate->Detector));
    }

    return s;
}

std::shared_ptr<ConstantRelativeEpow> ConstantRelativeEpow::make()
{
    return std::make_shared<ConstantRelativeEpow>();
}

std::shared_ptr<ConstantRelativeEpow> ConstantRelativeEpow::makeMCLongTarget()
{
    auto s = std::make_shared<ConstantRelativeEpow>();

    s->photon_cb   = { 0.0107, std_ext::degree_to_radian(3.79), std_ext::degree_to_radian(1.78)};
    s->photon_taps = { 0.035,  std_ext::degree_to_radian(0.42), std_ext::degree_to_radian(1.15)};

    s->proton_cb   = { 0.0,    std_ext::degree_to_radian(5.5), std_ext::degree_to_radian(5.3)};
    s->proton_taps = { 0.0,    std_ext::degree_to_radian(2.8), std_ext::degree_to_radian(4.45)};

    s->Eexp_cb   = -0.5;
    s->Eexp_taps = -0.5;

    return s;
}