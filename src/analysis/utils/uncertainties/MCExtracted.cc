#include "MCExtracted.h"

#include "expconfig/ExpConfig.h"

#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/math.h"
#include "base/std_ext/system.h"
#include "base/Paths.h"
#include "base/Logger.h"

#include "TH1D.h"
#include "TF1.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis::utils;
using namespace ant::analysis::utils::UncertaintyModels;

MCExtracted::angular_sigma::angular_sigma()
{}

MCExtracted::angular_sigma::~angular_sigma()
{}

double MCExtracted::angular_sigma::GetSigma(const unsigned element, const double E) const
{
    const double vp0 = p0->GetBinContent(int(element));
    const double vp1 = p1->GetBinContent(int(element));
    const double vp2 = p2->GetBinContent(int(element));

    return f(E, vp0, vp1, vp2);

}

double MCExtracted::angular_sigma::f(const double x, const double p0, const double p1, const double p2) noexcept
{
    //return p0*exp(p1*x)+p2;
    return exp(p0 + p1*x) + p2;
}

double MCExtracted::angular_sigma::f_root(const double* x, const double* p) noexcept
{
    return f(x[0], p[0], p[1], p[2]);
}

TF1* MCExtracted::angular_sigma::GetTF1(const std::string& name)
{
    auto f = new TF1(name.c_str(), f_root, 0, 1600, 3);
    f->SetParName(0, "#alpha");
    f->SetParName(1, "#beta");
    f->SetParName(2, "Offset");
    return f;
}




void MCExtracted::angular_sigma::Load(WrapTFile& f, const string& prefix, const int bins)
{
    p0 = LoadHist(f, prefix+"_p0", bins);
    p1 = LoadHist(f, prefix+"_p1", bins);
    p2 = LoadHist(f, prefix+"_p2", bins);
}

MCExtracted::angular_sigma::Hist MCExtracted::angular_sigma::LoadHist(WrapTFile& f, const string& name, const int bins)
{
    auto h = f.GetSharedHist<TH1D>(name);

    if(!h)
        throw std::runtime_error("Could not load histogram " + name);

    if(h->GetNbinsX() != bins)
        throw std::runtime_error(formatter() << "TH1D " << name << " has wrong number of bins: " << h->GetNbinsX()  << ", expected: " << bins);

    return h;
}

Uncertainties_t MCExtracted::GetSigmasProton(const TParticle &proton) const
{

    Uncertainties_t sigmas;

    sigmas.sigmaEk = 0.0;  // unmeasured

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

Uncertainties_t MCExtracted::GetSigmasPhoton(const TParticle &photon) const
{
    Uncertainties_t sigmas;

    const auto& cluster = photon.Candidate->FindCaloCluster();

    if(photon.Candidate->Detector & Detector_t::Type_t::CB) {

        sigmas.sigmaEk     = 1.07134e-02 * photon.Ek();
        sigmas.sigmaTheta = cb_sigma_theta.GetSigma(cluster->CentralElement, photon.Ek());
        sigmas.sigmaPhi   = cb_sigma_phi.GetSigma(cluster->CentralElement, photon.Ek());

    } else if(photon.Candidate->Detector & Detector_t::Type_t::TAPS) {

        sigmas.sigmaEk     = 3.5E-2 * photon.Ek();
        sigmas.sigmaTheta = taps_sigma_theta.GetSigma(cluster->CentralElement, photon.Ek());
        sigmas.sigmaPhi   = taps_sigma_phi.GetSigma(cluster->CentralElement, photon.Ek());

    } else {
        LOG(WARNING) << "Photon not in (CB,TAPS?)";
    }

    return sigmas;
}

MCExtracted::MCExtracted()
{}

MCExtracted::~MCExtracted()
{}

void MCExtracted::LoadSigmas(const string& filename)
{

    unique_ptr<WrapTFileInput> f;
    if(!std_ext::system::testopen(filename)) {
        LOG(WARNING) << "Could not read sigmas from '" << filename << "', using default.";
        f = std_ext::make_unique<WrapTFileInput>(string(ANT_PATH_DATABASE)+"/default/physics_files/FitterSigmas.root");
    }
    else {
        f = std_ext::make_unique<WrapTFileInput>(filename);
    }


    const auto& CB = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    const auto nCB = CB->GetNChannels();

    const auto& TAPS = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    const auto nTAPS = TAPS->GetNChannels();

    cb_sigma_theta.Load(*f, "CB_sigma_Theta", int(nCB));
    cb_sigma_phi.Load(  *f, "CB_sigma_Phi",   int(nCB));

    taps_sigma_theta.Load(*f, "TAPS_sigma_Theta", int(nTAPS));
    taps_sigma_phi.Load(  *f, "TAPS_sigma_Phi",   int(nTAPS));
}

Uncertainties_t MCExtracted::GetSigmas(const TParticle &particle) const
{

    if(particle.Type() == ParticleTypeDatabase::Photon)
        return GetSigmasPhoton(particle);

    if(particle.Type() == ParticleTypeDatabase::Proton)
        return GetSigmasProton(particle);

    // shoult never happen in normal running
    throw Exception("Unexpected Particle: " + particle.Type().Name());

}

std::shared_ptr<MCExtracted> MCExtracted::makeAndLoad()
{
    auto s = std::make_shared<MCExtracted>();
    auto& setup = ant::ExpConfig::Setup::Get();
    s->LoadSigmas(setup.GetPhysicsFilesDirectory()+"/FitterSigmas.root");
    return s;
}