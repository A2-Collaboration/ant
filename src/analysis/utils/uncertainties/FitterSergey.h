#pragma once

#include "analysis/utils/Uncertainties.h"

namespace ant {
namespace analysis {
namespace utils {
namespace UncertaintyModels {

/**
 * @brief Uncertainties from Sergey's fitter
 */
struct FitterSergey : public UncertaintyModel {
public:
    enum class beamtime_t {
        EPT_2014, Eta_2007
    };
    const beamtime_t Beamtime;

    FitterSergey(beamtime_t beamtime = beamtime_t::EPT_2014);
    virtual ~FitterSergey();

    Uncertainties_t GetSigmas(const TParticle& particle) const override;
};


}}}} // namespace ant::analysis::utils::UncertaintyModels