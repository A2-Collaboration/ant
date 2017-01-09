#pragma once

#include "analysis/utils/Uncertainties.h"

class TRandom;

namespace ant {
namespace analysis {
namespace utils {
namespace UncertaintyModels {

/**
 * @brief Uncertainties from Patrik Adlarson for MC Smearing
 *
 * It uses a random generator, for whatever reason.
 */
struct MCSmearingAdlarson : public UncertaintyModel {
public:

    MCSmearingAdlarson();
    virtual ~MCSmearingAdlarson();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<MCSmearingAdlarson> make();

protected:
    std::unique_ptr<TRandom> rng;
};

}}}} // namespace ant::analysis::utils::UncertaintyModels