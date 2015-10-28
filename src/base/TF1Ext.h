#pragma once

#include "TF1.h"

namespace ant {

template <typename T>
void FixParameters(TF1* func, const T& pars) {
    for(const auto& p : pars) {
        func->FixParameter(p, func->GetParameter(p));
    }
}

template <typename T>
void UnFixParameters(TF1* func, const T& pars) {
    for(const auto& p : pars) {
        func->ReleaseParameter(p);
    }
}

}
