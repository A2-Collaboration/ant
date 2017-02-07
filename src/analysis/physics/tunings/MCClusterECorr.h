#pragma once

#include "physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class MCClusterECorr : public Physics {

    const static BinSettings bins_Ek;

    TH2D* h_nCaloClusters = nullptr;
    TH2D* h_LostMCTrue = nullptr;

public:

    class CBTAPS_t {
        const Detector_t::Type_t Type;
        HistogramFactory HistFac;
        TH2D* h_nFills = nullptr;
        TH2D* h_EtrueErec = nullptr;
        TH3D* h_EtrueErec_3D = nullptr;
        TH2D* h_ErecEtrue_elements = nullptr;
    public:
        CBTAPS_t(Detector_t::Type_t type, const HistogramFactory& histFac,
                 const BinSettings& bins_cosTheta);
        void Fill(const TCluster& caloCluster, double Etrue) const;
        void Finish() const;
        void Draw(canvas& c) const;
    };

private:

    const CBTAPS_t CB;
    const CBTAPS_t TAPS;

public:
    MCClusterECorr(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};


}}} // end of ant::analysis::physics