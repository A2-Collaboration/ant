#include "Optimized.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/TAPS.h"
#include "base/std_ext/math.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "cereal/archives/json.hpp"
#pragma GCC diagnostic pop

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis::utils;
using namespace ant::analysis::utils::UncertaintyModels;

Optimized::Optimized()
{}

Optimized::~Optimized()
{}

Uncertainties_t Optimized::GetSigmas(const TParticle& particle) const
{
    if(!particle.Candidate)
        throw Exception("No candidate attached to particle");

    auto calocluster = particle.Candidate->FindCaloCluster();

    if(!calocluster)
        throw Exception("No calo cluster found");

    const auto theta = particle.Theta();
    const auto Ek     = particle.Ek();

    Uncertainties_t s;

    if(particle.Candidate->Detector & Detector_t::Type_t::CB) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            s.sigmaEk     = dE(Ek, cb_photon_E_rel, cb_photon_E_exp, cb_photon_E_lin);
            s.sigmaTheta = dThetaSin(theta, cb_photon_theta_const, cb_photon_theta_Sin);
            s.sigmaPhi   = cb_photon_phi / sin(theta);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            s = cb_proton;
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

        static auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
        auto elem = cb->GetClusterElement(calocluster->CentralElement);
        s.ShowerDepth = elem->RadiationLength*std::log2(Ek/elem->CriticalE); // /std::pow(std::sin(theta),3.0);
        s.sigmaCB_R = 2; // in cm
    }
    else if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            s.sigmaEk     = dE(Ek, taps_photon_E_rel, taps_photon_E_exp, taps_photon_E_lin);
            s.sigmaTheta = taps_photon_theta;
            s.sigmaPhi   = taps_photon_phi;
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            s = taps_proton;
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

        static auto taps = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
        auto elem = taps->GetClusterElement(calocluster->CentralElement);
        s.ShowerDepth = elem->RadiationLength*std::log2(Ek/elem->CriticalE);
        s.sigmaTAPS_L  = 2; // in cm
        s.sigmaTAPS_Rxy = 1; // in cm
    }
    else {
        throw Exception("Unexpected Detector: " + string(particle.Candidate->Detector));
    }

    return s;
}

double Optimized::dThetaSin(const double theta, const double offset, const double thetapart) noexcept
{
    return offset + thetapart * sin(theta);
}

double Optimized::dE(const double E, const double rel, const double exp, const double reloffset) noexcept
{
    return rel * E * pow( E/1000.0, exp) + reloffset * E;
}

string angleoutput(const double x) {
    const auto y = radian_to_degree(x);
    return formatter() << setprecision(3) << y;
}

string numberoutput(const double x) {
    return formatter() << setprecision(4) << x;
}

double angleinput(const string& x) {
    return degree_to_radian(atof(x.c_str()));
}

double numberinput(const string& x) {
    return atof(x.c_str());
}

string Optimized::to_string_simple() const
{
    return formatter()
            << "cgtc="<< angleoutput(cb_photon_theta_const)  << separator
            << "cgts="<< angleoutput(cb_photon_theta_Sin)    << separator
            << "cgp=" << angleoutput(cb_photon_phi)          << separator
            << "cgEr="<< numberoutput(cb_photon_E_rel)       << separator
            << "cgEe="<< numberoutput(cb_photon_E_exp)       << separator
            << "cgEl="<< numberoutput(cb_photon_E_lin)       << separator
            << "cpt=" << angleoutput(cb_proton.sigmaTheta)   << separator
            << "cpp=" << angleoutput(cb_proton.sigmaPhi)     << separator
            << "tgt=" << angleoutput(taps_photon_theta)      << separator
            << "tgp=" << angleoutput(taps_photon_phi)        << separator
            << "tgEr="<< numberoutput(taps_photon_E_rel)     << separator
            << "tgEe="<< numberoutput(taps_photon_E_exp)     << separator
            << "tgEl="<< numberoutput(taps_photon_E_lin)     << separator
            << "tpt=" << angleoutput(taps_proton.sigmaTheta) << separator
            << "tpp=" << angleoutput(taps_proton.sigmaPhi);
}

string Optimized::to_string() const
{
    std::stringstream ss; // any stream can be used

    {
      cereal::JSONOutputArchive oarchive(ss); // Create an output archive

      const auto cb_photon_theta_const_d = angleoutput(cb_photon_theta_const);
      const auto cb_photon_theta_Sin_d   = angleoutput(cb_photon_theta_Sin);
      const auto cb_photon_phi_d         = angleoutput(cb_photon_phi);
      const auto cb_photon_E_rel_d       = numberoutput(cb_photon_E_rel);
      const auto cb_photon_E_exp_d       = numberoutput(cb_photon_E_exp);
      const auto cb_photon_E_lin_d       = numberoutput(cb_photon_E_lin);
      const auto cb_proton_theta_d       = angleoutput(cb_proton.sigmaTheta);
      const auto cb_proton_phi_d         = angleoutput(cb_proton.sigmaPhi);
      const auto taps_photon_theta_d     = angleoutput(taps_photon_theta);
      const auto taps_photon_phi_d       = angleoutput(taps_photon_phi);
      const auto taps_photon_E_rel_d     = numberoutput(taps_photon_E_rel);
      const auto taps_photon_E_exp_d     = numberoutput(taps_photon_E_exp);
      const auto taps_photon_E_lin_d     = numberoutput(taps_photon_E_lin);
      const auto taps_proton_theta_d     = angleoutput(taps_proton.sigmaTheta);
      const auto taps_proton_phi_d       = angleoutput(taps_proton.sigmaPhi);


      oarchive(
                  cereal::make_nvp("cgtc", cb_photon_theta_const_d),
                  cereal::make_nvp("cgts", cb_photon_theta_Sin_d),
                  cereal::make_nvp("cgp",  cb_photon_phi_d),
                  cereal::make_nvp("cgEr", cb_photon_E_rel_d),
                  cereal::make_nvp("cgEe", cb_photon_E_exp_d),
                  cereal::make_nvp("cgEl", cb_photon_E_lin_d),
                  cereal::make_nvp("cpt",  cb_proton_theta_d),
                  cereal::make_nvp("cpp",  cb_proton_phi_d),
                  cereal::make_nvp("tgt",  taps_photon_theta_d),
                  cereal::make_nvp("tgp",  taps_photon_phi_d),
                  cereal::make_nvp("tgEr", taps_photon_E_rel_d),
                  cereal::make_nvp("tgEe", taps_photon_E_exp_d),
                  cereal::make_nvp("tgEl", taps_photon_E_lin_d),
                  cereal::make_nvp("tpt",  taps_proton_theta_d),
                  cereal::make_nvp("tpp",  taps_proton_phi_d)
                  );
    }

    return ss.str();
}

string Optimized::to_string_short() const
{
    auto str = to_string();
    str.erase(std::remove_if(str.begin(),
                                  str.end(),
                                  [](char x){return std::isspace(x);}),
                   str.end());
    return str;
}

void Optimized::load_from_string(const string& data)
{
    std::stringstream ss(data);

    cereal::JSONInputArchive iarchive(ss);

    string cb_photon_theta_const_d;
    string cb_photon_theta_Sin_d;
    string cb_photon_phi_d;
    string cb_photon_E_rel_d;
    string cb_photon_E_exp_d;
    string cb_photon_E_lin_d;
    string cb_proton_theta_d;
    string cb_proton_phi_d;
    string taps_photon_theta_d;
    string taps_photon_phi_d;
    string taps_photon_E_rel_d;
    string taps_photon_E_exp_d;
    string taps_photon_E_lin_d;
    string taps_proton_theta_d;
    string taps_proton_phi_d;

    iarchive(
                cereal::make_nvp("cgtc", cb_photon_theta_const_d),
                cereal::make_nvp("cgts", cb_photon_theta_Sin_d),
                cereal::make_nvp("cgp",  cb_photon_phi_d),
                cereal::make_nvp("cgEr", cb_photon_E_rel_d),
                cereal::make_nvp("cgEe", cb_photon_E_exp_d),
                cereal::make_nvp("cgEe", cb_photon_E_lin_d),
                cereal::make_nvp("cpt",  cb_proton_theta_d),
                cereal::make_nvp("cpp",  cb_proton_phi_d),
                cereal::make_nvp("tgt",  taps_photon_theta_d),
                cereal::make_nvp("tgp",  taps_photon_phi_d),
                cereal::make_nvp("tgEr", taps_photon_E_rel_d),
                cereal::make_nvp("tgEe", taps_photon_E_exp_d),
                cereal::make_nvp("tgEl", taps_photon_E_lin_d),
                cereal::make_nvp("tpt",  taps_proton_theta_d),
                cereal::make_nvp("tpp",  taps_proton_phi_d)
                );

    cb_photon_theta_const  = angleinput(cb_photon_theta_const_d);
    cb_photon_theta_Sin    = angleinput(cb_photon_theta_Sin_d);
    cb_photon_phi          = angleinput(cb_photon_phi_d);
    cb_photon_E_rel        = numberinput(cb_photon_E_rel_d);
    cb_photon_E_exp        = numberinput(cb_photon_E_exp_d);
    cb_photon_E_lin        = numberinput(cb_photon_E_lin_d);
    cb_proton.sigmaEk       = 0.0;
    cb_proton.sigmaTheta   = angleinput(cb_proton_theta_d);
    cb_proton.sigmaPhi     = angleinput(cb_proton_phi_d);
    taps_photon_theta      = angleinput(taps_photon_theta_d);
    taps_photon_phi        = angleinput(taps_photon_phi_d);
    taps_photon_E_rel      = numberinput(taps_photon_E_rel_d);
    taps_photon_E_exp      = numberinput(taps_photon_E_exp_d);
    taps_photon_E_lin      = numberinput(taps_photon_E_lin_d);
    taps_proton.sigmaEk     = 0.0;
    taps_proton.sigmaTheta = angleinput(taps_proton_theta_d);
    taps_proton.sigmaPhi   = angleinput(taps_proton_phi_d);
}

void Optimized::load_from_string_simple(const string& data)
{

    const auto tokens = std_ext::tokenize_string(data, separator);

    for(const auto& token : tokens) {
        ReadToken(token);
    }

}


void Optimized::ReadToken(const string& token)
{
    const auto tokens = std_ext::tokenize_string(token, "=");
    if(tokens.size() != 2)
        return;

    const auto& name  = tokens.front();
    const auto& value = tokens.back();

    if(name == "cgtc") {
        cb_photon_theta_const = angleinput(value);
    } else if(name=="cgts") {
        cb_photon_theta_Sin = angleinput(value);
    } else if(name=="cgp") {
        cb_photon_phi = angleinput(value);
    } else if(name=="cgEr") {
        cb_photon_E_rel = numberinput(value);
    } else if(name=="cgEl") {
        cb_photon_E_lin = numberinput(value);
    } else if(name=="cgEe") {
        cb_photon_E_exp = numberinput(value);
    } else if(name=="cpt") {
        cb_proton.sigmaTheta = angleinput(value);
    } else if(name=="cpp") {
        cb_proton.sigmaPhi = angleinput(value);
    } else if(name=="tgt") {
        taps_photon_theta = angleinput(value);
    } else if(name=="tgp") {
        taps_photon_phi = angleinput(value);
    } else if(name=="tgEr") {
        taps_photon_E_rel = numberinput(value);
    } else if(name=="tgEe") {
        taps_photon_E_exp = numberinput(value);
    } else if(name=="tgEl") {
        taps_photon_E_lin = numberinput(value);
    } else if(name=="tpt") {
        taps_proton.sigmaTheta = angleinput(value);
    } else if(name=="tpp") {
        taps_proton.sigmaPhi = angleinput(value);
    }
}

Optimized_Oli1::Optimized_Oli1(double relative_scale, bool use_measured_proton_TAPS)
{
    cb_photon_theta_const = relative_scale*degree_to_radian(1.1);
    cb_photon_theta_Sin   = relative_scale*degree_to_radian(3.9);
    cb_photon_phi         = relative_scale*degree_to_radian(4.1); // as average from dTheta over theta

    cb_photon_E_rel       =  relative_scale*0.02;  // 2% of E
    cb_photon_E_exp       = -0.25;  // dev by 4th root of E (in GeV)
    cb_photon_E_lin       =  0.0;   // no linear part

    cb_proton   = { 0.0, relative_scale*degree_to_radian(5.5), relative_scale*degree_to_radian(5.3)};

    taps_photon_E_rel =  relative_scale*0.03;    // 3% of E
    taps_photon_E_exp = -0.5;     // dev by sqrt
    taps_photon_E_lin =  relative_scale*0.018;   // 1.8% of E as linear part

    taps_photon_theta = relative_scale*degree_to_radian(2.5);
    taps_photon_phi   = relative_scale*degree_to_radian(2.0);


    taps_proton = { use_measured_proton_TAPS ? 30.0 : 0.0, relative_scale*degree_to_radian(2.8), relative_scale*degree_to_radian(4.45)};

}

Optimized_Oli1::~Optimized_Oli1()
{

}

const std::string Optimized::separator = ": ";