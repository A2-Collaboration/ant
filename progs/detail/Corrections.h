#pragma once

#include <vector>
#include <algorithm>

#include "base/vec/LorentzVec.h"

#include "TLorentzVector.h"
#include "PParticle.h"

#include "dilepton_radiative_corrections.h"
#include "dimuon_radiative_corrections.h"

namespace ant {
namespace progs {
namespace corrections {

enum class Channel { pi0_eeg, eta_eeg, eta_mumug, etap_eeg, etap_mumug, unknown };

double interpolate_linear(const double value, const double x_low, const double x_up, const double y_low, const double y_up);
double interpolate_bilinear(const double x, const double y, const double x_low, const double x_up, const double y_low, const double y_up,
                            const double val_x_low_y_low, const double val_x_up_y_low, const double val_x_low_y_up, const double val_x_up_y_up);

// helper class to apply radiative corrections to Dalitz decays
class Corrections {
  public:
    Corrections() = delete;
    Corrections(const Channel);
    ~Corrections() = default;
    Corrections(const Corrections&) = default;
    Corrections& operator=(const Corrections&) = delete;

    double interpolate(const double, const double) const;
    double interpolate(const TLorentzVector* const, const TLorentzVector* const, const TLorentzVector* const) const;
    double interpolate(const PParticle* const, const PParticle* const, const PParticle* const) const;

    double get_max_correction() const;

  private:
    const Channel channel;
    // vectors containing data points to determine the limits of the correction table
    std::vector<double>* y_vals = nullptr;
    // pointer to the x and y values as well as the corrections
    // used as a fallback to approximate the correction if the delauny interpolation fails
    std::vector<double>* x_vec = nullptr;
    std::vector<double>* y_vec = nullptr;
    std::vector<double>* z_vec = nullptr;
};

}}} // namespace ant::progs::corrections
