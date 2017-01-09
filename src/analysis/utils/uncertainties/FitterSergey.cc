#include "FitterSergey.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"


using namespace std;
using namespace ant;
using namespace ant::analysis::utils;
using namespace ant::analysis::utils::UncertaintyModels;

UncertaintyModels::FitterSergey::FitterSergey(beamtime_t beamtime) :
    Beamtime(beamtime)
{

}

UncertaintyModels::FitterSergey::~FitterSergey()
{

}

Uncertainties_t UncertaintyModels::FitterSergey::GetSigmas(const TParticle& particle) const
{
    if(!particle.Candidate)
        throw Exception("No candidate attached to particle");

    auto calocluster = particle.Candidate->FindCaloCluster();

    if(!calocluster)
        throw Exception("No calo cluster found");

    const auto& Theta = particle.Theta();
    const auto& Ek     = particle.Ek();

    Uncertainties_t u;
    u.Detector = particle.Candidate->Detector;

    constexpr double GeVMeV = 1000.0;

    if(u.Detector & Detector_t::Type_t::CB)
    {
        if(particle.Type() == ParticleTypeDatabase::Photon) {

            auto dEovEclCB = [this] (double Ecl) {
                auto dEovEclCBInit = [] (double Ecl) {
                    const double p[5] = {0.014, 0.0025, 0.35, 0., 0.0032};
                    return p[0] / pow(Ecl + p[1], p[2]) + p[3] + p[4] * Ecl;
                };
                auto dEovEclCBAdd = [this] (double Ecl) {
                    (void)Ecl; // prevent unused variable warning
                    switch(Beamtime) {
                    case beamtime_t::EPT_2014:
                        return 0.052;
                    case beamtime_t::Eta_2007:
                        return 0.0255;
                    }
                    return std_ext::NaN;
                };
                double Er0 = dEovEclCBInit(Ecl);
                double Era = dEovEclCBAdd(Ecl);
                return sqrt(pow(Er0, 2) + pow(Era, 2));
            };

            auto dThetaCB = [] (double Ecl) {
                const double p[4] = {7.69518e-03, 4.86197e-01, 1.79483, 1.57948e-02};
                return p[0] / pow(Ecl + p[1], p[2]) + p[3];
            };

            auto DepthShowCB = [] (double Ecl) {
                const double p[4] = {-3.36631, 9.40334e-02, 5.35372e-01, 4.36397e+01};
                return p[0] / pow(Ecl + p[1], p[2]) + p[3];
            };

            auto dDepthShowCB = [] (double Ecl) {
                const double p[4] = {1.76634e-01, 0., 6.26983e-01, 2.48218};
                double sdep = p[0] / pow(Ecl + p[1], p[2]) + p[3];
                sdep *= 1.05;
                return sdep;
            };

            u.sigmaEk = dEovEclCB(Ek/GeVMeV)*Ek;
            u.sigmaTheta = dThetaCB(Ek/GeVMeV);
            u.ShowerDepth = DepthShowCB(Ek/GeVMeV);
            u.sigmaCB_R = dDepthShowCB(Ek/GeVMeV);
        }
        else if(particle.Type() == ParticleTypeDatabase::Proton) {
            auto dThetaCB = [] (double Ecl) {
                double p[4] = {1.38476e-04, 5.30098e-01, 7.61558, 3.75841e-02};
                p[3] += 0.004;
                double dTh = p[0] / pow(Ecl + p[1], p[2]) + p[3];
                dTh *= 1.25;
                return dTh;
            };

            auto DepthShowCB = [] (double Ecl) {
               const double p[4] = {2.52512e+01, 6.44248, 1.96292e+02, -1.61958e+02};
               return p[0] + Ecl * p[1] + Ecl * Ecl * p[2] + Ecl * Ecl * Ecl * p[3];
            };

            auto dDepthShowCB = [] (double Ecl) {
                const double p[4] = {3.5783e-02, 3.47172e-01, 1.50307, -4.88434e-01};
                double sdep = p[0] + Ecl * p[1] + Ecl * Ecl * p[2] + Ecl * Ecl * Ecl * p[3];
                sdep *= 1.05;
                return sdep;
            };

            u.sigmaEk = 0; // proton Ek is unmeasured
            u.sigmaTheta = dThetaCB(Ek/GeVMeV);
            u.ShowerDepth = DepthShowCB(Ek/GeVMeV);
            u.sigmaCB_R = dDepthShowCB(Ek/GeVMeV);

        }
        else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

        // Sergey's CB shower depths include the inner CB radius, but Ant uncertainties do not
        static auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
        u.ShowerDepth -= cb->GetInnerRadius();

        u.sigmaPhi = u.sigmaTheta / sin(Theta);
    }
    else if(u.Detector & Detector_t::Type_t::TAPS)
    {
        if(particle.Type() == ParticleTypeDatabase::Photon) {

            auto dEovEclTAPS = [] (double Ecl) {
                auto dEovEclTAPSInit = [] (double Ecl) {
                    double p[4] = {1.88319e-04, 1.42657, 3.96356e-02, 1.52351e-02};
                    p[3] *= 1.8;
                    double Er0 = p[0] / pow(Ecl - 0.002, p[1]) + p[2] + p[3] * Ecl;
                    return Er0;
                };
                auto dEovEclTAPSAdd = [] (double Ecl) {
                    return 0.031 + 0.04 * Ecl;
                };
                double Er0 = dEovEclTAPSInit(Ecl);
                double Era = dEovEclTAPSAdd(Ecl);
                return sqrt(pow(Er0, 2) + pow(Era, 2));
            };

            auto dTanThTAPS = [] (double Ecl) {
                const double p[5] = {3.28138e+02, 0., 7.29002e-04, -3.27381e+02, 0.};
                double dtan = p[0] / pow(Ecl + p[1], p[2]) + p[3];
                dtan *= 0.85;
                return dtan;
            };

            auto DepthShowTAPS = [] (double Ecl) {
                double p[4] = {-2.99791e+01, 1.75852e-03, 4.99643e-02, 4.14362e+01};
                p[0] *= 0.978;
                return p[0] / pow(Ecl + p[1], p[2]) + p[3];
            };

            auto dDepthShowTAPS = [] (double Ecl) {
                const double p[4] = {2.83139, 0., 1.02537e-01, -7.53507e-01};
                return p[0] / pow(Ecl + p[1], p[2]) + p[3];
            };

            u.sigmaEk = dEovEclTAPS(Ek/GeVMeV)*Ek;
            u.sigmaTAPS_Rxy = dTanThTAPS(Ek/GeVMeV);
            u.ShowerDepth = DepthShowTAPS(Ek/GeVMeV);
            u.sigmaTAPS_L = dDepthShowTAPS(Ek/GeVMeV);

            auto& pos = particle.Candidate->FindCaloCluster()->Position;
            const auto& TAPS_Lz = pos.z + u.ShowerDepth*std::cos(particle.Theta());
            const vec3 TAPS_L{pos.x, pos.y, TAPS_Lz};
            const auto& TAPS_Rxy = std::sin(particle.Theta())*TAPS_L.R();

            u.sigmaPhi = u.sigmaTAPS_Rxy / TAPS_Rxy;
        }
        else if(particle.Type() == ParticleTypeDatabase::Proton) {

            auto dTanThTAPS = [] (double Ecl, double RadCl) {
                const double p[4] = {3.27709e+02, 4.99670e-02, 5.55520e-03, -3.27819e+02};
                double dtan = p[0] / pow(Ecl + p[1], p[2]) + p[3];
                if (RadCl > 41.)
                    dtan *= 1.3;
                return dtan;
            };

            auto DepthShowTAPS = [] (double Ecl) {
                const double p[4] = {-1.73216e-02, 3.83753, 1.54891e+02, -1.328e+02};
                double dprot = p[0] + Ecl * p[1] + Ecl * Ecl * p[2] + Ecl * Ecl * Ecl * p[3];
                return dprot * 1.05;
            };

            auto dDepthShowTAPS = [] (double Ecl) {
                const double p[4] = {8.43187e-03, 3.63264e-01, 7.17476e-01,  7.33715};
                return p[0] + Ecl * p[1] + Ecl * Ecl * p[2] + Ecl * Ecl * Ecl * p[3];
            };

            u.ShowerDepth = DepthShowTAPS(Ek/GeVMeV);

            auto& pos = particle.Candidate->FindCaloCluster()->Position;
            const auto& TAPS_Lz = pos.z + u.ShowerDepth*std::cos(particle.Theta());
            const vec3 TAPS_L{pos.x, pos.y, TAPS_Lz};
            const auto& TAPS_Rxy = std::sin(particle.Theta())*TAPS_L.R();

            u.sigmaEk = 0; // proton is unmeasured
            u.sigmaTAPS_Rxy = dTanThTAPS(Ek/GeVMeV, TAPS_Rxy);
            u.sigmaTAPS_L = dDepthShowTAPS(Ek/GeVMeV);
            u.sigmaPhi = u.sigmaTAPS_Rxy / TAPS_Rxy;
        }
        else {
            throw Exception("Unexpected Particle: " + particle.Type().Name());
        }

    }
    else {
        throw Exception("Unexpected Detector: " + string(particle.Candidate->Detector));
    }


    return u;
}

double UncertaintyModels::FitterSergey::GetBeamEnergySigma(double) const
{
    return 1.3;
}