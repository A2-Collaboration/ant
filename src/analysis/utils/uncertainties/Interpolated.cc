#include "Interpolated.h"

#include "expconfig/ExpConfig.h"

#include "base/Interpolator.h"
#include "base/Array2D.h"
#include "base/std_ext/memory.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/ClippedInterpolatorWrapper.h"

#include "TAxis.h"
#include "TH2D.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;
using namespace ant::analysis::utils::UncertaintyModels;

Interpolated::Interpolated(UncertaintyModelPtr starting_uncertainty_, bool use_proton_sigmaE_) :
    starting_uncertainty(starting_uncertainty_),
    use_proton_sigmaE(use_proton_sigmaE_)
{
    if(use_proton_sigmaE && !starting_uncertainty)
        LOG(WARNING) << "Requested to use proton sigmaE from starting uncertainty model, but no such model provided";
}

Interpolated::~Interpolated()
{

}

Uncertainties_t Interpolated::GetSigmas(const TParticle& particle) const
{
    auto u_starting = starting_uncertainty ?
                          starting_uncertainty->GetSigmas(particle)
                        : Uncertainties_t{};

    if(!loaded_sigmas)
        return u_starting;

    auto u = u_starting;
    auto& detector = particle.Candidate->Detector;
    if(detector & Detector_t::Type_t::CB) {
        if(particle.Type() == ParticleTypeDatabase::Photon) {
            cb_photon.SetUncertainties(u, particle);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            cb_proton.SetUncertainties(u, particle);
        } else {
            throw Exception("Unexpected Particle in CB: " + particle.Type().Name());
        }
    } else if(detector & Detector_t::Type_t::TAPS) {
        if(particle.Type() == ParticleTypeDatabase::Photon) {
            taps_photon.SetUncertainties(u, particle);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            taps_proton.SetUncertainties(u, particle);
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }
    }
    else {
        throw Exception("Unexpected Detector: " + string(detector));
    }

    // special handling for proton E uncertainty (is unmeasured if not flagged)
    if(particle.Type() == ParticleTypeDatabase::Proton) {
        //
        if(starting_uncertainty && use_proton_sigmaE) {
            u.sigmaEk = u_starting.sigmaEk;
        }
        // sanitize interpolation of "zero" uncertainty in any case
        if(!std::isfinite(u.sigmaEk) || u.sigmaEk < 1e-5) {
            u.sigmaEk = 0;
        }
    }

    return u;
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


void Interpolated::LoadSigmas(const string& filename)
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
        LOG(WARNING) << "Can't open uncertainty histogram file (using default model instead): " << e.what();
    }

}



std::shared_ptr<Interpolated> Interpolated::makeAndLoad(
        Type_t type,
        UncertaintyModelPtr default_model,
        bool use_proton_sigmaE)
{
    auto s = std::make_shared<Interpolated>(default_model, use_proton_sigmaE);

    auto& setup = ant::ExpConfig::Setup::Get();

    s->LoadSigmas(
                std_ext::formatter()
                << setup.GetPhysicsFilesDirectory()
                << "/interpolated_sigmas" << (type == Type_t::MC ? "_mc" : "_data") << ".root"
                );
    if(!default_model && !s->HasLoadedSigmas()) {
        throw Exception("No default model provided and sigmas could not be loaded");
    }
    return s;
}

namespace ant {
namespace analysis {
namespace utils {
namespace UncertaintyModels {

ostream& operator<<(ostream& stream, const Interpolated& o)
{
    stream << "Photon CB:\n"   << o.cb_photon   << "\n";
    stream << "Proton CB:\n"   << o.cb_proton   << "\n";
    stream << "Photon TAPS:\n" << o.taps_photon << "\n";
    stream << "Proton TAPS:\n" << o.taps_proton << "\n";
    return stream;
}

ostream& operator<<(ostream& stream, const Interpolated::EkThetaPhiR& o)
{
    stream << "Ek:\t\t"        << o.Ek     << "\n";
    stream << "Theta:\t\t"     << o.Theta << "\n";
    stream << "Phi:\t\t"       << o.Phi   << "\n";
    stream << "CB_R:\t\t"      << o.CB_R  << "\n";
    stream << "ShowerDepth:\t" << o.ShowerDepth  << "\n";
    return stream;
}

ostream& operator<<(ostream& stream, const Interpolated::EkRxyPhiL& o)
{
    stream << "Ek:\t\t"        << o.Ek     << "\n";
    stream << "TAPS_Rxy:\t"    << o.TAPS_Rxy << "\n";
    stream << "Phi:\t\t"       << o.Phi   << "\n";
    stream << "TAPS_L:\t\t"    << o.TAPS_L  << "\n";
    stream << "ShowerDepth:\t" << o.ShowerDepth  << "\n";
    return stream;
}

}}}} // namespace ant::analysis::utils::UncertaintyModels

std::unique_ptr<const Interpolator2D> Interpolated::LoadInterpolator(const WrapTFile& file, const string& hname)
{
    TH2D* h = nullptr;

    file.GetObject(hname, h);

    if(!h)
        throw(std::runtime_error("Histogram not found: "+hname));

    return ClippedInterpolatorWrapper::makeInterpolator(h);

}

void Interpolated::EkThetaPhiR::SetUncertainties(Uncertainties_t& u, const TParticle& particle) const
{
    auto costheta = std::cos(particle.Theta());
    auto Ekin = particle.Ek();
    u.sigmaEk     = Ek.GetPoint(costheta, Ekin);
    u.sigmaTheta  = Theta.GetPoint(costheta, Ekin);
    u.sigmaPhi    = Phi.GetPoint(costheta, Ekin);
    u.sigmaCB_R   = CB_R.GetPoint(costheta, Ekin);
    u.ShowerDepth = ShowerDepth.GetPoint(costheta, Ekin);
}

void Interpolated::EkThetaPhiR::Load(const WrapTFile& file, const std::string& prefix)
{
    Ek.setInterpolator(    LoadInterpolator(file, prefix+"/sigma_Ek"));
    Theta.setInterpolator( LoadInterpolator(file, prefix+"/sigma_Theta"));
    Phi.setInterpolator(   LoadInterpolator(file, prefix+"/sigma_Phi"));
    CB_R.setInterpolator(  LoadInterpolator(file, prefix+"/sigma_R"));
    ShowerDepth.setInterpolator(  LoadInterpolator(file, prefix+"/h_NewShowerDepth"));
}

void Interpolated::EkRxyPhiL::SetUncertainties(Uncertainties_t& u, const TParticle& particle) const
{
    auto costheta = std::cos(particle.Theta());
    auto Ekin = particle.Ek();
    u.sigmaEk       = Ek.GetPoint(costheta, Ekin);
    u.sigmaTAPS_Rxy = TAPS_Rxy.GetPoint(costheta, Ekin);
    u.sigmaPhi      = Phi.GetPoint(costheta, Ekin);
    u.sigmaTAPS_L   = TAPS_L.GetPoint(costheta, Ekin);
    u.ShowerDepth   = ShowerDepth.GetPoint(costheta, Ekin);
}

void Interpolated::EkRxyPhiL::Load(const WrapTFile& file, const std::string& prefix)
{
    Ek.setInterpolator(       LoadInterpolator(file, prefix+"/sigma_Ek"));
    TAPS_Rxy.setInterpolator( LoadInterpolator(file, prefix+"/sigma_Rxy"));
    Phi.setInterpolator(      LoadInterpolator(file, prefix+"/sigma_Phi"));
    TAPS_L.setInterpolator(   LoadInterpolator(file, prefix+"/sigma_L"));
    ShowerDepth.setInterpolator(   LoadInterpolator(file, prefix+"/h_NewShowerDepth"));
}

