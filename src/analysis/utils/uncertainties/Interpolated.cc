#include "Interpolated.h"

#include "expconfig/ExpConfig.h"

#include "base/Interpolator.h"
#include "base/Array2D.h"
#include "base/std_ext/memory.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"

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

}

Interpolated::~Interpolated()
{

}

Uncertainties_t Interpolated::GetSigmas(const TParticle& particle) const
{
    auto u = starting_uncertainty ? starting_uncertainty->GetSigmas(particle) : Uncertainties_t{};

    if(!loaded_sigmas)
        return u;

    if(u.Detector & Detector_t::Type_t::CB) {
        if(particle.Type() == ParticleTypeDatabase::Photon) {
            cb_photon.SetUncertainties(u, particle);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            cb_proton.SetUncertainties(u, particle);
        } else {
            throw Exception("Unexpected Particle in CB: " + particle.Type().Name());
        }
    } else if(u.Detector & Detector_t::Type_t::TAPS) {
        if(particle.Type() == ParticleTypeDatabase::Photon) {
            taps_photon.SetUncertainties(u, particle);
        } else if(particle.Type() == ParticleTypeDatabase::Proton) {
            taps_proton.SetUncertainties(u, particle);
        } else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }
    }
    else {
        throw Exception("Unexpected Detector: " + string(particle.Candidate->Detector));
    }

    // special handling for proton E uncertainty (is unmeasured if not flagged)
    if(particle.Type() == ParticleTypeDatabase::Proton) {
        if (
           !use_proton_sigmaE
           || !std::isfinite(u.sigmaEk) || u.sigmaEk < 1e-5  // sanitize interpolation
           )
        {
            if(starting_uncertainty) {
                // fallback to starting_model, but only for energy!
                auto u_starting = starting_uncertainty->GetSigmas(particle);
                u.sigmaEk = u_starting.sigmaEk;
            }
            else {
                u.sigmaEk = 0;
            }
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

struct OrderedGrid2D {
    vector<double> x;
    vector<double> y;
    Array2D        z;

    OrderedGrid2D(const unsigned nx, const unsigned ny, const double dflt=0.0):
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

    OrderedGrid2D grid(nx+pad_x*2, ny+pad_y*2, std_ext::NaN);

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

    // flood fill averages the rest
    FloodFillAverages::fillNeighborAverages(grid.z);

    return std_ext::make_unique<Interpolator2D>(grid.x,grid.y, grid.z.Data());
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

    return makeInterpolator(h);

}

void Interpolated::EkThetaPhiR::SetUncertainties(Uncertainties_t& u, const TParticle& particle) const
{
    auto costheta = std::cos(particle.Theta());
    auto Ekin = particle.Ek();
    u.sigmaEk    = Ek.GetPoint(costheta, Ekin);
    u.sigmaTheta = Theta.GetPoint(costheta, Ekin);
    u.sigmaPhi   = Phi.GetPoint(costheta, Ekin);
    u.sigmaCB_R  = CB_R.GetPoint(costheta, Ekin);
}

void Interpolated::EkThetaPhiR::Load(const WrapTFile& file, const std::string& prefix)
{
    Ek.setInterpolator(    LoadInterpolator(file, prefix+"/sigma_Ek"));
    Theta.setInterpolator( LoadInterpolator(file, prefix+"/sigma_Theta"));
    Phi.setInterpolator(   LoadInterpolator(file, prefix+"/sigma_Phi"));
    CB_R.setInterpolator(  LoadInterpolator(file, prefix+"/sigma_R"));
}

ostream& Interpolated::EkThetaPhiR::Print(ostream& stream) const
{
    stream << "Ek:\t\t"    << Ek     << "\n";
    stream << "Theta:\t\t" << Theta << "\n";
    stream << "Phi:\t\t"   << Phi   << "\n";
    stream << "CB_R:\t\t"  << CB_R  << "\n";
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
}

void Interpolated::EkRxyPhiL::Load(const WrapTFile& file, const std::string& prefix)
{
    Ek.setInterpolator(       LoadInterpolator(file, prefix+"/sigma_Ek"));
    TAPS_Rxy.setInterpolator( LoadInterpolator(file, prefix+"/sigma_Rxy"));
    Phi.setInterpolator(      LoadInterpolator(file, prefix+"/sigma_Phi"));
    TAPS_L.setInterpolator(   LoadInterpolator(file, prefix+"/sigma_L"));
}

ostream& Interpolated::EkRxyPhiL::Print(ostream& stream) const
{
    stream << "Ek:\t\t"       << Ek     << "\n";
    stream << "TAPS_Rxy:\t"   << TAPS_Rxy << "\n";
    stream << "Phi:\t\t"      << Phi   << "\n";
    stream << "TAPS_L:\t\t"   << TAPS_L  << "\n";
    return stream;
}


Interpolated::ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(Interpolated::ClippedInterpolatorWrapper::interpolator_ptr_t i):
    interp(move(i)), xrange(interp->getXRange()), yrange(interp->getYRange())
{}

Interpolated::ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(): interp(nullptr), xrange({0,0}), yrange({0,0})
{}

Interpolated::ClippedInterpolatorWrapper::~ClippedInterpolatorWrapper()
{}

double Interpolated::ClippedInterpolatorWrapper::GetPoint(double x, double y) const
{
    x = xrange.clip(x);
    y = yrange.clip(y);
    return interp->GetPoint(x,y);
}

void Interpolated::ClippedInterpolatorWrapper::setInterpolator(Interpolated::ClippedInterpolatorWrapper::interpolator_ptr_t i) {
    interp = move(i);
    xrange = interp->getXRange();
    yrange = interp->getYRange();
}

ostream& Interpolated::ClippedInterpolatorWrapper::Print(ostream& stream) const
{
    stream << "Theta Range: " << xrange << " E Range: " << yrange;
    return stream;
}


double Interpolated::ClippedInterpolatorWrapper::boundsCheck_t::clip(double v) const
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

ostream& Interpolated::ClippedInterpolatorWrapper::boundsCheck_t::Print(ostream& stream) const
{
    stream << range << "-> [" << underflow << "|" << unclipped << "|" << overflow << "]";
    return stream;
}