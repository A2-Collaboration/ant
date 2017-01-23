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
        UncertaintyModelPtr default_model,
        bool use_proton_sigmaE)
{
    auto s = std::make_shared<Interpolated>(default_model, use_proton_sigmaE);

    const auto setup = ant::ExpConfig::Setup::GetLastFound();

    if(!setup) {
        throw ExpConfig::ExceptionNoConfig("No Setup found");
    }

    s->LoadSigmas(
                std_ext::formatter()
                << setup->GetPhysicsFilesDirectory()
                << "/interpolated_sigmas.root"
                );
    if(!default_model && !s->HasLoadedSigmas()) {
        throw Exception("No default model provided and sigmas could not be loaded");
    }
    return s;
}

ostream& Interpolated::Print(ostream& stream) const
{
    stream << "Photon CB:\n"   << cb_photon   << "\n";
    stream << "Proton CB:\n"   << cb_proton   << "\n";
    stream << "Photon TAPS:\n" << taps_photon << "\n";
    stream << "Proton TAPS:\n" << taps_proton << "\n";
    return stream;
}

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

ostream& Interpolated::EkThetaPhiR::Print(ostream& stream) const
{
    stream << "Ek:\t\t"        << Ek     << "\n";
    stream << "Theta:\t\t"     << Theta << "\n";
    stream << "Phi:\t\t"       << Phi   << "\n";
    stream << "CB_R:\t\t"      << CB_R  << "\n";
    stream << "ShowerDepth:\t" << ShowerDepth  << "\n";
    return stream;
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

ostream& Interpolated::EkRxyPhiL::Print(ostream& stream) const
{
    stream << "Ek:\t\t"        << Ek     << "\n";
    stream << "TAPS_Rxy:\t"    << TAPS_Rxy << "\n";
    stream << "Phi:\t\t"       << Phi   << "\n";
    stream << "TAPS_L:\t\t"    << TAPS_L  << "\n";
    stream << "ShowerDepth:\t" << ShowerDepth  << "\n";
    return stream;
}
