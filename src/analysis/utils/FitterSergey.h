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

    struct result_t : printable_traits {
        double   TaggE;
        double   TaggT;
        unsigned TaggCh;

        double KinFitProb;
        double TreeFitProb;
        double AntiPi0FitProb;
        double AntiEtaFitProb;

        std::vector<double> IM_3g;
        double IM_4g;

        std::vector<double> gNonPi0_Theta;
        std::vector<double> gNonPi0_CaloE;

        double CBVetoSumE;

        virtual std::ostream& Print(std::ostream& stream) const;
    };

    std::vector<result_t> Process(const TEventData& data);

};



}}} // namespace ant::analysis::utils