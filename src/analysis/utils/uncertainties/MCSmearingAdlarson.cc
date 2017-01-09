#include "MCSmearingAdlarson.h"

#include "base/std_ext/memory.h"

#include "TRandom2.h"
#include "TMath.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;
using namespace ant::analysis::utils::UncertaintyModels;

MCSmearingAdlarson::MCSmearingAdlarson() :
    rng(std_ext::make_unique<TRandom2>())
{

}

MCSmearingAdlarson::~MCSmearingAdlarson()
{

}

Uncertainties_t MCSmearingAdlarson::GetSigmas(const TParticle& particle) const
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

std::shared_ptr<MCSmearingAdlarson> MCSmearingAdlarson::make()
{
    return std::make_shared<MCSmearingAdlarson>();
}