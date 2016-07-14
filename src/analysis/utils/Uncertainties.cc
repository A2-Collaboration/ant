#include "Uncertainties.h"

#include "expconfig/ExpConfig.h"

#include "base/WrapTFile.h"
#include "base/Paths.h"
#include "base/Logger.h"
#include "base/std_ext/system.h"
#include "base/Interpolator.h"
#include "base/Array2D.h"

#include <ostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "base/cereal/archives/json.hpp"
#pragma GCC diagnostic pop

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom2.h"
#include "TMath.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;

UncertaintyModel::~UncertaintyModel()
{}


UncertaintyModels::MCExtracted::angular_sigma::angular_sigma()
{}

UncertaintyModels::MCExtracted::angular_sigma::~angular_sigma()
{}

double UncertaintyModels::MCExtracted::angular_sigma::GetSigma(const unsigned element, const double E) const
{
    const double vp0 = p0->GetBinContent(int(element));
    const double vp1 = p1->GetBinContent(int(element));
    const double vp2 = p2->GetBinContent(int(element));

    return f(E, vp0, vp1, vp2);

}

double UncertaintyModels::MCExtracted::angular_sigma::f(const double x, const double p0, const double p1, const double p2) noexcept
{
    //return p0*exp(p1*x)+p2;
    return exp(p0 + p1*x) + p2;
}

double UncertaintyModels::MCExtracted::angular_sigma::f_root(const double* x, const double* p) noexcept
{
    return f(x[0], p[0], p[1], p[2]);
}

TF1* UncertaintyModels::MCExtracted::angular_sigma::GetTF1(const std::string& name)
{
    auto f = new TF1(name.c_str(), f_root, 0, 1600, 3);
    f->SetParName(0, "#alpha");
    f->SetParName(1, "#beta");
    f->SetParName(2, "Offset");
    return f;
}




void UncertaintyModels::MCExtracted::angular_sigma::Load(WrapTFile& f, const string& prefix, const int bins)
{
    p0 = LoadHist(f, prefix+"_p0", bins);
    p1 = LoadHist(f, prefix+"_p1", bins);
    p2 = LoadHist(f, prefix+"_p2", bins);
}

UncertaintyModels::MCExtracted::angular_sigma::Hist UncertaintyModels::MCExtracted::angular_sigma::LoadHist(WrapTFile& f, const string& name, const int bins)
{
    auto h = f.GetSharedHist<TH1D>(name);

    if(!h)
        throw std::runtime_error("Could not load histogram " + name);

    if(h->GetNbinsX() != bins)
        throw std::runtime_error(formatter() << "TH1D " << name << " has wrong number of bins: " << h->GetNbinsX()  << ", expected: " << bins);

    return h;
}

Uncertainties_t UncertaintyModels::MCExtracted::GetSigmasProton(const TParticle &proton) const
{

    Uncertainties_t sigmas;

    sigmas.sigmaE = 0.0;  // unmeasured

    if(proton.Candidate->Detector & Detector_t::Type_t::CB) {
        sigmas.sigmaTheta = degree_to_radian(5.43);
        sigmas.sigmaPhi   = degree_to_radian(5.31);
    } else if(proton.Candidate->Detector & Detector_t::Type_t::TAPS) {
        sigmas.sigmaTheta = degree_to_radian(2.89);
        sigmas.sigmaPhi   = degree_to_radian(4.45);
    } else {
        LOG(WARNING) << "Proton is not in (CB,TAPS)?";
    }

    return sigmas;
}

Uncertainties_t UncertaintyModels::MCExtracted::GetSigmasPhoton(const TParticle &photon) const
{
    Uncertainties_t sigmas;

    const auto& cluster = photon.Candidate->FindCaloCluster();

    if(photon.Candidate->Detector & Detector_t::Type_t::CB) {

        sigmas.sigmaE     = 1.07134e-02 * photon.Ek();
        sigmas.sigmaTheta = cb_sigma_theta.GetSigma(cluster->CentralElement, photon.Ek());
        sigmas.sigmaPhi   = cb_sigma_phi.GetSigma(cluster->CentralElement, photon.Ek());

    } else if(photon.Candidate->Detector & Detector_t::Type_t::TAPS) {

        sigmas.sigmaE     = 3.5E-2 * photon.Ek();
        sigmas.sigmaTheta = taps_sigma_theta.GetSigma(cluster->CentralElement, photon.Ek());
        sigmas.sigmaPhi   = taps_sigma_phi.GetSigma(cluster->CentralElement, photon.Ek());

    } else {
        LOG(WARNING) << "Photon not in (CB,TAPS?)";
    }

    return sigmas;
}

UncertaintyModels::MCExtracted::MCExtracted()
{}

UncertaintyModels::MCExtracted::~MCExtracted()
{}

void UncertaintyModels::MCExtracted::LoadSigmas(const string& filename)
{
    const auto& setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw Exception("No Setup found!");

    unique_ptr<WrapTFileInput> f;
    if(!std_ext::system::testopen(filename)) {
        LOG(WARNING) << "Could not read sigmas from '" << filename << "', using default.";
        f = std_ext::make_unique<WrapTFileInput>(string(ANT_PATH_DATABASE)+"/default/physics_files/FitterSigmas.root");
    }
    else {
        f = std_ext::make_unique<WrapTFileInput>(filename);
    }


    const auto& CB = setup->GetDetector(Detector_t::Type_t::CB);
    if(!CB)
        throw Exception("No CB detector defined in setup");

    const auto nCB = CB->GetNChannels();

    const auto& TAPS = setup->GetDetector(Detector_t::Type_t::TAPS);
    if(!TAPS)
        throw Exception("No TAPS detector defined in setup");

    const auto nTAPS = TAPS->GetNChannels();

    cb_sigma_theta.Load(*f, "CB_sigma_Theta", int(nCB));
    cb_sigma_phi.Load(  *f, "CB_sigma_Phi",   int(nCB));

    taps_sigma_theta.Load(*f, "TAPS_sigma_Theta", int(nTAPS));
    taps_sigma_phi.Load(  *f, "TAPS_sigma_Phi",   int(nTAPS));
}

Uncertainties_t UncertaintyModels::MCExtracted::GetSigmas(const TParticle &particle) const
{

    if(particle.Type() == ParticleTypeDatabase::Photon)
        return GetSigmasPhoton(particle);

    if(particle.Type() == ParticleTypeDatabase::Proton)
        return GetSigmasProton(particle);

    // shoult never happen in normal running
    throw Exception("Unexpected Particle: " + particle.Type().Name());

}

std::shared_ptr<UncertaintyModels::MCExtracted> UncertaintyModels::MCExtracted::makeAndLoad()
{
    auto s = std::make_shared<MCExtracted>();

    const auto setup = ant::ExpConfig::Setup::GetLastFound();

    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

    s->LoadSigmas(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    return s;
}


UncertaintyModels::Constant::Constant()
{}

UncertaintyModels::Constant::~Constant()
{}

std::shared_ptr<UncertaintyModels::Constant> UncertaintyModels::Constant::make()
{
    return std::make_shared<Constant>();
}

Uncertainties_t ant::analysis::utils::UncertaintyModels::Constant::GetSigmas(const TParticle &particle) const
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

UncertaintyModels::ConstantRelativeE::ConstantRelativeE()
{}

UncertaintyModels::ConstantRelativeE::~ConstantRelativeE()
{}

Uncertainties_t UncertaintyModels::ConstantRelativeE::GetSigmas(const TParticle &particle) const
{
    auto s = Constant::GetSigmas(particle);
    s.sigmaE *= particle.Ek();

    return s;
}

std::shared_ptr<UncertaintyModels::ConstantRelativeE> UncertaintyModels::ConstantRelativeE::makeMCLongTarget()
{
    auto s = std::make_shared<ConstantRelativeE>();

    s->photon_cb   = { 0.0107, std_ext::degree_to_radian(3.79), std_ext::degree_to_radian(1.78)};
    s->photon_taps = { 0.035,  std_ext::degree_to_radian(0.42), std_ext::degree_to_radian(1.15)};

    s->proton_cb   = { 0.0,    std_ext::degree_to_radian(5.5), std_ext::degree_to_radian(5.3)};
    s->proton_taps = { 0.0,    std_ext::degree_to_radian(2.8), std_ext::degree_to_radian(4.45)};

    return s;
}

UncertaintyModels::ConstantRelativeEpow::ConstantRelativeEpow()
{}

UncertaintyModels::ConstantRelativeEpow::~ConstantRelativeEpow()
{}

Uncertainties_t UncertaintyModels::ConstantRelativeEpow::GetSigmas(const TParticle &particle) const
{
    Uncertainties_t s;

    if(particle.Candidate->Detector & Detector_t::Type_t::CB) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            s = photon_cb;
            s.sigmaE = s.sigmaE* particle.Ek() * pow(particle.Ek(), Eexp_cb);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            s = proton_cb;
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

    } else if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {
            s = photon_taps;
            s.sigmaE = s.sigmaE* particle.Ek() * pow(particle.Ek(), Eexp_taps);
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

std::shared_ptr<UncertaintyModels::ConstantRelativeEpow> UncertaintyModels::ConstantRelativeEpow::make()
{
    return std::make_shared<ConstantRelativeEpow>();
}

std::shared_ptr<UncertaintyModels::ConstantRelativeEpow> UncertaintyModels::ConstantRelativeEpow::makeMCLongTarget()
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

UncertaintyModels::Optimized::Optimized()
{}

UncertaintyModels::Optimized::~Optimized()
{}

Uncertainties_t UncertaintyModels::Optimized::GetSigmas(const TParticle& particle) const
{
    const auto theta = particle.Theta();
    const auto E     = particle.Ek();

    Uncertainties_t s;

    if(particle.Candidate->Detector & Detector_t::Type_t::CB) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {

            s.sigmaE     = dE(E, cb_photon_E_rel, cb_photon_E_exp, cb_photon_E_lin);
            s.sigmaTheta = dThetaSin(theta, cb_photon_theta_const, cb_photon_theta_Sin);
            s.sigmaPhi   = cb_photon_phi / sin(theta);

        } else if(particle.Type() == ParticleTypeDatabase::Proton) {

            s = cb_proton;

        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

    } else if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {

            s.sigmaE     = dE(E, taps_photon_E_rel, taps_photon_E_exp, taps_photon_E_lin);
            s.sigmaTheta = taps_photon_theta;
            s.sigmaPhi   = taps_photon_phi;

        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            s = taps_proton;

        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }
    }
    else {
        throw Exception("Unexpected Detector: " + string(particle.Candidate->Detector));
    }

    return s;
}

double UncertaintyModels::Optimized::dThetaSin(const double theta, const double offset, const double thetapart) noexcept
{
    return offset + thetapart * sin(theta);
}

double UncertaintyModels::Optimized::dE(const double E, const double rel, const double exp, const double reloffset) noexcept
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

string UncertaintyModels::Optimized::to_string_simple() const
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

string UncertaintyModels::Optimized::to_string() const
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

string UncertaintyModels::Optimized::to_string_short() const
{
    auto str = to_string();
    str.erase(std::remove_if(str.begin(),
                                  str.end(),
                                  [](char x){return std::isspace(x);}),
                   str.end());
    return str;
}

void UncertaintyModels::Optimized::load_from_string(const string& data)
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
    cb_proton.sigmaE       = 0.0;
    cb_proton.sigmaTheta   = angleinput(cb_proton_theta_d);
    cb_proton.sigmaPhi     = angleinput(cb_proton_phi_d);
    taps_photon_theta      = angleinput(taps_photon_theta_d);
    taps_photon_phi        = angleinput(taps_photon_phi_d);
    taps_photon_E_rel      = numberinput(taps_photon_E_rel_d);
    taps_photon_E_exp      = numberinput(taps_photon_E_exp_d);
    taps_photon_E_lin      = numberinput(taps_photon_E_lin_d);
    taps_proton.sigmaE     = 0.0;
    taps_proton.sigmaTheta = angleinput(taps_proton_theta_d);
    taps_proton.sigmaPhi   = angleinput(taps_proton_phi_d);
}

void UncertaintyModels::Optimized::load_from_string_simple(const string& data)
{

    const auto tokens = std_ext::tokenize_string(data, separator);

    for(const auto& token : tokens) {
        ReadToken(token);
    }

}

bool UncertaintyModels::Optimized::operator==(const UncertaintyModels::Optimized& other) const noexcept
{
    return
               cb_photon_theta_const == other.cb_photon_theta_const
            && cb_photon_theta_Sin   == other.cb_photon_theta_Sin
            && cb_photon_phi         == other.cb_photon_phi
            && cb_photon_E_rel       == other.cb_photon_E_rel
            && cb_photon_E_exp       == other.cb_photon_E_exp
            && cb_photon_E_lin       == other.cb_photon_E_lin
            && cb_proton             == other.cb_proton
            && taps_photon_E_rel     == other.taps_photon_E_rel
            && taps_photon_E_exp     == other.taps_photon_E_exp
            && taps_photon_E_lin     == other.taps_photon_E_lin
            && taps_photon_theta     == other.taps_photon_theta
            && taps_photon_phi       == other.taps_photon_phi
            && taps_proton           == other.taps_proton;

}

bool UncertaintyModels::Optimized::operator!=(const UncertaintyModels::Optimized& other) const noexcept
{
    return !(*this == other);
}

void UncertaintyModels::Optimized::ReadToken(const string& token)
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

UncertaintyModels::Optimized_Oli1::Optimized_Oli1()
{
    cb_photon_theta_const = degree_to_radian(1.1);
    cb_photon_theta_Sin   = degree_to_radian(3.9);
    cb_photon_phi         = degree_to_radian(4.1); // as average from dTheta over theta

    cb_photon_E_rel       =  0.02;  // 2% of E
    cb_photon_E_exp       = -0.25;  // dev by 4th root of E (in GeV)
    cb_photon_E_lin       =  0.0;   // no linear part

    cb_proton   = { 0.0, degree_to_radian(5.5), degree_to_radian(5.3)};

    taps_photon_E_rel =  0.03;    // 3% of E
    taps_photon_E_exp = -0.5;     // dev by sqrt
    taps_photon_E_lin =  0.018;   // 1.8% of E as linear part

    taps_photon_theta = degree_to_radian(2.5);
    taps_photon_phi   = degree_to_radian(2.0);


    taps_proton = { 0.0, degree_to_radian(2.8), degree_to_radian(4.45)};

}

const std::string UncertaintyModels::Optimized::separator = ": ";


UncertaintyModels::MCSmearingAdlarson::MCSmearingAdlarson() :
    rng(std_ext::make_unique<TRandom2>())
{

}

UncertaintyModels::MCSmearingAdlarson::~MCSmearingAdlarson()
{

}

Uncertainties_t UncertaintyModels::MCSmearingAdlarson::GetSigmas(const TParticle& particle) const
{
    if(particle.Candidate->Detector & Detector_t::Type_t::CB) {
        // code taken from
        // https://github.com/padlarson/acqu/blob/work/acqu_user/root/src/TA2CalArray.cc#L233
        // parameters taken from
        // https://github.com/padlarson/acqu/blob/work/acqu_user/data.MC/AR-Analysis-CentApp-NaI.dat#L32

        constexpr double MC_Smear_ThetaMin = 20.0;
        constexpr double MC_Smear_ThetaMax = 160.0;
        constexpr double MC_Smear_CBThetaBoundary = 28.0;
        constexpr double MC_Smear_ThetaSigma = 2.5;
        constexpr double MC_SmearMin = 0.010;
        constexpr double MC_SmearMax = 0.060;
        constexpr double MC_Smear_EnergyMax = 1.2;
        constexpr double MC_SmearBoundaryMin = 0.001;
        constexpr double MC_SmearBoundaryMax = 0.1;



        double theta = std_ext::radian_to_degree(particle.Theta());
        auto energy = particle.Ek() / 1000; // kinetic energy in GeV


        // smear theta angle
        theta += rng->Gaus(0.0, MC_Smear_ThetaSigma);

        // "decay" constant to mimic the experimental resolution
        double c = ( TMath::Log(MC_SmearMax/MC_SmearMin) )/( MC_Smear_ThetaMax-MC_Smear_ThetaMin );
        // calculate smearing value, mutliply by the minimum smearing value
        double sigma = MC_SmearMin*TMath::Exp(c*(MC_Smear_ThetaMax-theta));

        // increase smearing closer to the CB boundary region
        if( theta < MC_Smear_CBThetaBoundary && energy < MC_Smear_EnergyMax){
            // linearly decreasing in E
            sigma += MC_SmearBoundaryMax - energy*(MC_SmearBoundaryMax - MC_SmearBoundaryMin)/MC_Smear_EnergyMax;
        }

        // https://github.com/padlarson/acqu/blob/work/acqu_user/data.MC/AR-Analysis-CentApp-NaI.dat#L20
        // from D.Werthmueller PhD thesis...
        constexpr double fSigmaEnergyPower = 0.66;

        sigma *= TMath::Power(energy, fSigmaEnergyPower);

        return {sigma*1000, 0, 0}; // energy smearing only, convert back to MeV
    }
    else if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {
        // code taken from
        // https://github.com/padlarson/acqu/blob/work/acqu_user/root/src/TA2TAPS_BaF2.cc#L278
        // parameters taken from
        // https://github.com/padlarson/acqu/blob/work/acqu_user/data.MC/AR-Analysis-TAPS-BaF2.dat#L31

        constexpr double MC_SmearMin = 0.03;
        constexpr double MC_SmearMax = 0.07;
        constexpr double MC_Smear_EnergyMax = 0.5;
        constexpr double MC_Smear_EnergyRes = 0.01;

        double sigma; // smearing to be applied, in GeV
        auto energy = particle.Ek() / 1000; // kinetic energy in GeV

        if(energy < MC_Smear_EnergyMax){
            sigma = MC_SmearMax - energy/MC_Smear_EnergyMax * (MC_SmearMax - MC_SmearMin);
        }
        else
            sigma = 0;

        sigma = sigma*sqrt(energy) + MC_Smear_EnergyRes*energy;

        return {sigma*1000, 0, 0}; // energy smearing only, convert back to MeV
    }

    // no smearing by default
    return {0,0,0};
}

std::shared_ptr<UncertaintyModels::MCSmearingAdlarson> UncertaintyModels::MCSmearingAdlarson::make()
{
    return std::make_shared<UncertaintyModels::MCSmearingAdlarson>();
}




UncertaintyModels::Interpolated::Interpolated(UncertaintyModelPtr starting_uncertainty_) :
    starting_uncertainty(starting_uncertainty_)
{

}

UncertaintyModels::Interpolated::~Interpolated()
{

}

Uncertainties_t UncertaintyModels::Interpolated::GetSigmas(const TParticle& particle) const
{
    // use starting model if nothing was loaded so far
    if(!loaded_sigmas) {
        return starting_uncertainty->GetSigmas(particle);
    }


    if(particle.Candidate->Detector & Detector_t::Type_t::CB) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {

            return cb_photon.GetUncertainties(particle);

        } else if(particle.Type() == ParticleTypeDatabase::Proton) {

            auto u = cb_proton.GetUncertainties(particle);
            u.sigmaE = 0.0;
            return u;

        } else {
            throw Exception("Unexpected Particle in CB: " + particle.Type().Name());
        }

    } else if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {

        if(particle.Type() == ParticleTypeDatabase::Photon) {

            return taps_photon.GetUncertainties(particle);

        } else if(particle.Type() == ParticleTypeDatabase::Proton) {

            auto u = taps_proton.GetUncertainties(particle);
            u.sigmaE = 0.0;
            return u;

        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }
    }
    else {
        throw Exception("Unexpected Detector: " + string(particle.Candidate->Detector));
    }
}

std::vector<double> getBinPositions(const TAxis* axis) {

    assert(axis->GetNbins() >= 1);

    vector<double> bins(size_t(axis->GetNbins())+2);

    for(int i=1; i<=axis->GetNbins(); ++i) {
        bins.at(size_t(i)) = axis->GetBinCenter(i);
    }

    const auto binwidth = axis->GetBinWidth(1); //assume equal bin sizes everywhere

    bins.front() = bins.at(1) - binwidth;
    bins.back()  = bins.at(bins.size()-2) + binwidth;

    return bins;
}

struct OrnderedGrid2D {
    vector<double> x;
    vector<double> y;
    Array2D        z;

    OrnderedGrid2D(const unsigned nx, const unsigned ny, const double dflt=0.0):
        x(nx),
        y(ny),
        z(nx, ny, dflt) {}

};

std::vector<double> getBinContents(const TH2D* hist) {
    const auto nx = hist->GetNbinsX();
    const auto ny = hist->GetNbinsY();

    assert(nx>0);
    assert(ny>0);

    vector<double> z(size_t(nx*ny));

    for(int x=1; x<=nx;++x) {
        for(int y=1; y<ny+1; ++y) {
            z.at(size_t((x-1)+nx*(y-1))) = hist->GetBinContent(x,y);
        }
    }

    return z;
}

std::unique_ptr<const Interpolator2D> makeInterpolator(TH2D* hist) {

    const unsigned nx = unsigned(hist->GetNbinsX());
    const unsigned pad_x = nx > 3 ? 1 : 2;

    const unsigned ny = unsigned(hist->GetNbinsY());
    const unsigned pad_y = ny > 3 ? 1 : 2;

    OrnderedGrid2D grid(nx+pad_x*2, ny+pad_y*2, std_ext::NaN);

    // extend x bin positions
    {
        double avgBinWidth = .0;
        for(unsigned x=0; x<nx; ++x) {
            grid.x.at(x+pad_x) = hist->GetXaxis()->GetBinCenter(int(x+1));
            avgBinWidth += hist->GetXaxis()->GetBinWidth(int(x+1));
        }
        avgBinWidth /= nx;
        for(unsigned x=1; x<=pad_x;++x) {
            grid.x.at(pad_x-x)    = grid.x.at(pad_x-x+1)    - avgBinWidth;
            grid.x.at(pad_x+nx+x-1) = grid.x.at(pad_x+nx+x-2) + avgBinWidth;
        }
    }

    // eytend y bin positions
    {
        double avgBinWidth = .0;
        for(unsigned y=0; y<ny; ++y) {
            grid.y.at(y+pad_y) = hist->GetYaxis()->GetBinCenter(int(y+1));
            avgBinWidth += hist->GetYaxis()->GetBinWidth(int(y+1));
        }
        avgBinWidth /= ny;
        for(unsigned y=1; y<=pad_y;++y) {
            grid.y.at(pad_y-y)    = grid.y.at(pad_y-y+1)    - avgBinWidth;
            grid.y.at(pad_y+ny+y-1) = grid.y.at(pad_y+ny+y-2) + avgBinWidth;
        }
    }

    // copy data to the "middle". leave some borders around
    grid.z.CopyRect(Array2D_TH2D(hist), pad_x, pad_y);



    // fill borders:

    // top and botton rows
    for(unsigned y=1; y<=pad_y; ++y ) {
        for(unsigned x=0; x<nx; ++x) {
            grid.z.at(x+pad_x,pad_y-y)       = grid.z.at(x+pad_x,pad_y-y+1);
            grid.z.at(x+pad_x,pad_y+ny+y-1)    = grid.z.at(x+pad_x,pad_y+ny+y-2);
        }
    }

    // first and last column
    for(unsigned x=1; x<=pad_x; ++x ) {
        for(unsigned y=0; y<ny; ++y) {
            grid.z.at(pad_x-x,   y+pad_y)      = grid.z.at(pad_x-x+1,   y+pad_y);
            grid.z.at(pad_x+nx+x-1,y+pad_y)    = grid.z.at(pad_x+nx+x-2,y+pad_y);
        }
    }

    //flood fill averages the rest
    FloodFillAverages::fillNeighborAverages(grid.z);



    return std_ext::make_unique<Interpolator2D>(grid.x,grid.y, grid.z.Data());
}

void UncertaintyModels::Interpolated::LoadSigmas(const string& filename)
{
    try {

        WrapTFileInput f(filename);

        cb_photon.Load  (f, "sigma_photon_cb");
        cb_proton.Load  (f, "sigma_proton_cb");
        taps_photon.Load(f, "sigma_photon_taps");
        taps_proton.Load(f, "sigma_proton_taps");

        loaded_sigmas = true;
        VLOG(5) << "Successfully loaded interpolation data for Uncertainty Model from " << filename;

    } catch (WrapTFile::Exception& e) {
        LOG(WARNING) << "Can't open uncertainty histogram file: " << filename <<". Default model will be used: " << e.what();
    }

}



std::shared_ptr<UncertaintyModels::Interpolated> UncertaintyModels::Interpolated::makeAndLoad(UncertaintyModelPtr default_model)
{
    auto s = std::make_shared<Interpolated>(default_model);

    const auto setup = ant::ExpConfig::Setup::GetLastFound();

    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

    s->LoadSigmas(setup->GetPhysicsFilesDirectory()+"/interpolated_sigmas.root");

    return s;
}

ostream&UncertaintyModels::Interpolated::Print(ostream& stream) const
{
    stream << "Photon CB:\n"   << cb_photon   << "\n";
    stream << "Proton CB:\n"   << cb_proton   << "\n";
    stream << "Photon TAPS:\n" << taps_photon << "\n";
    stream << "Proton TAPS:\n" << taps_proton << "\n";
    return stream;
}

Uncertainties_t UncertaintyModels::Interpolated::EkThetaPhi::GetUncertainties(const TParticle& particle) const
{
    auto costheta = std::cos(particle.Theta());
    auto Ek = particle.Ek();

    return {
        E.GetPoint(costheta, Ek),
        Theta.GetPoint(costheta, Ek),
        Phi.GetPoint(costheta, Ek)
    };
}

void UncertaintyModels::Interpolated::EkThetaPhi::Load(WrapTFile& file, const std::string& prefix)
{
    E.setInterpolator(LoadInterpolator(file, prefix+"/sigma_E"));
    Theta.setInterpolator(LoadInterpolator(file, prefix+"/sigma_Theta"));
    Phi.setInterpolator(LoadInterpolator(file, prefix+"/sigma_Phi"));

}

std::unique_ptr<const Interpolator2D> UncertaintyModels::Interpolated::EkThetaPhi::LoadInterpolator(WrapTFile& file, const string& hname)
{
    TH2D* h = nullptr;

    file.GetObject(hname, h);

    if(!h)
        throw(std::runtime_error("Histogram not found: "+hname));

    return makeInterpolator(h);

}

ostream&UncertaintyModels::Interpolated::EkThetaPhi::Print(ostream& stream) const
{
    stream << "E:\t"     << E     << "\n";
    stream << "Theta:\t" << Theta << "\n";
    stream << "Phi:\t"   << Phi   << "\n";
    return stream;
}

UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::interpolator_ptr_t i):
    interp(move(i)), xrange(interp->getXRange()), yrange(interp->getYRange())
{}

UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(): interp(nullptr), xrange({0,0}), yrange({0,0})
{}

UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::~ClippedInterpolatorWrapper()
{}

double UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::GetPoint(double x, double y) const
{
    x = xrange.clip(x);
    y = yrange.clip(y);
    return interp->GetPoint(x,y);
}

void UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::setInterpolator(UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::interpolator_ptr_t i) {
    interp = move(i);
    xrange = interp->getXRange();
    yrange = interp->getYRange();
}

ostream&UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::Print(ostream& stream) const
{
    stream << "Theta Range: " << xrange << " E Range: " << yrange;
    return stream;
}


double UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::boundsCheck_t::clip(double v) const
{
    if(v < range.Start()) {
        underflow++;
        return range.Start();
    }

    if(v > range.Stop()) {
        overflow++;
        return range.Stop();
    }

    unclipped++;

    return v;
}

ostream& UncertaintyModels::Interpolated::ClippedInterpolatorWrapper::boundsCheck_t::Print(ostream& stream) const
{
    stream << range << "-> [" << underflow << "|" << unclipped << "|" << overflow << "]";
    return stream;
}



