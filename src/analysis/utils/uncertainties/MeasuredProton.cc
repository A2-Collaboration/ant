#include "MeasuredProton.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;
using namespace ant::analysis::utils::UncertaintyModels;

MeasuredProton::MeasuredProton(std::shared_ptr<const UncertaintyModel>& base_):
    UncertaintyModel(),
    base(base_)
{

}

MeasuredProton::~MeasuredProton()
{

}

Uncertainties_t MeasuredProton::GetSigmas(const TParticle& particle) const
{

    auto res = base->GetSigmas(particle);

    if(particle.Type() == ParticleTypeDatabase::Proton) {

        if(!particle.Candidate)
            throw Exception("No candidate attached to particle");

        const auto& Ek     = particle.Ek();

        constexpr double GeVMeV = 1000.0;


        auto dEovEclCBInit = [] (double Ecl) { // CB dE/E

            double p[5] = {5.69464e-05, 1.48943e-01, 3.41725, 1.11244e-02,
                             -1.77329e-03};

            p[0] = 0.043;
            p[1] = 0.;
            p[2] = 0.43;
            p[3] = 0.;
            p[4] = 0.;

            return p[0] / pow(Ecl + p[1], p[2]) + p[3] + p[4] * Ecl;
        };

        auto dEovEclTAPSInit = [] (double Ecl) { // TAPS dE/E

            double p[4] = {1.88319e-04, 1.42657, 3.96356e-02, 1.52351e-02};

            p[0] = 0.045;
            p[1] = 0.;
            p[2] = 0.45;
            p[3] = 0.;
            return p[0] / pow(Ecl + p[1], p[2]) + p[3];
        };

        if(particle.Candidate->Detector & Detector_t::Type_t::CB) {
            res.sigmaEk = dEovEclCBInit(Ek/GeVMeV)*Ek;
        } else  if(particle.Candidate->Detector & Detector_t::Type_t::TAPS) {
            res.sigmaEk = dEovEclTAPSInit(Ek/GeVMeV)*Ek;
        }

    }

    return res;
}