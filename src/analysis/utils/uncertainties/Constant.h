#pragma once

#include "analysis/utils/Uncertainties.h"

namespace ant {
namespace analysis {
namespace utils {
namespace UncertaintyModels {

/**
 * @brief Constant/static kin fitter uncertainty model. has a fixed value for {cb,taps}{photon,proton}{E,theta,phi}
 */
struct Constant : public UncertaintyModel {
public:

    Uncertainties_t photon_cb;
    Uncertainties_t photon_taps;
    Uncertainties_t proton_cb;
    Uncertainties_t proton_taps;

    Constant();
    virtual ~Constant();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<Constant> make();
};

/**
 * @brief Simple kin fitter uncertainty model. has a fixed value for {cb,taps}{photon,proton}{theta,phi},
 * Energies are relative values and get multipied with the particle energy on GetSigmas()
 */
struct ConstantRelativeE : public Constant {
public:

    ConstantRelativeE();
    virtual ~ConstantRelativeE();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new, empty instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeE> make();

    /**
     * @brief Create a new instance filled with global values determined from MC and return a shared ptr to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeE> makeMCLongTarget();
};

struct ConstantRelativeEpow : public Constant {
public:
    double Eexp_cb   = 1.0;
    double Eexp_taps = 1.0;

    ConstantRelativeEpow();
    virtual ~ConstantRelativeEpow();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new, empty instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeEpow> make();

    /**
     * @brief Create a new instance filled with global values determined from MC and return a shared ptr to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeEpow> makeMCLongTarget();
};

}}}} // namespace ant::analysis::utils::UncertaintyModels