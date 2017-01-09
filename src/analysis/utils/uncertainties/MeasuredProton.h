#pragma once

#include "analysis/utils/Uncertainties.h"

namespace ant {
namespace analysis {
namespace utils {
namespace UncertaintyModels {

struct MeasuredProton : public UncertaintyModel {
    std::shared_ptr<const UncertaintyModel> base = nullptr;

    MeasuredProton(std::shared_ptr<const UncertaintyModel>& base_);
    virtual ~MeasuredProton();

    Uncertainties_t GetSigmas(const TParticle& particle) const override;
};

}}}} // namespace ant::analysis::utils::UncertaintyModels