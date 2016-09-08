#pragma once

#include "tree/TEventData.h"

#include <memory>

namespace ant {
namespace analysis {
namespace utils {

class FitterSergey {

    // use PIMPL idiom to hide all the implementation (copied from Sergey's Acqu)
    class TA2KFitC;
    std::unique_ptr<TA2KFitC> I;

public:
    FitterSergey();
    virtual ~FitterSergey();

    struct result_t {
        double TaggT;
        double TaggCh;

    };

    std::vector<result_t> Process(const TEventData& data);

};



}}} // namespace ant::analysis::utils