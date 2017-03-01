#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/root_draw.h"

#include "analysis/utils/A2GeoAcceptance.h"

class TH3;
class TH1D;

namespace ant {
namespace analysis {
namespace physics {


class GeoAcceptance: public Physics {

private:

    class ParticleThetaPhiPlot: public ant::root_drawable_traits {
    public:
        TH2D* hist;
        const static BinSettings theta_bins;
        const static BinSettings phi_bins;

        ParticleThetaPhiPlot(HistogramFactory& factory,
                             const std::string& title,
                             const std::string& name="",
                             const BinSettings& thetabins = BinSettings(180),
                             const BinSettings& phibins = BinSettings(360,-180.0,180.0));
        void Fill(const TParticlePtr& p);


        TObject *GetObject();
        void Draw(const std::string &option) const;
    };

    class ParticleThetaPhiPlot3D: public ant::root_drawable_traits {
    public:
        TH3* hist;
        unsigned int n;

        ParticleThetaPhiPlot3D(HistogramFactory& factory, const std::string& title, const std::string& name="");
        void Fill(const TParticlePtr& p);


        TObject *GetObject();
        void Draw(const std::string &option) const;
    };

    class AcceptanceAnalysis {
    public:
        std::string name;
        HistogramFactory HistFac;
        const utils::A2SimpleGeometry& geo;
        ParticleThetaPhiPlot mctrue_pos;
        ParticleThetaPhiPlot matched_pos;
        ParticleThetaPhiPlot multimatched_pos;
        ParticleThetaPhiPlot lost_pos;
        ParticleThetaPhiPlot lost_pos_zoom;
        ParticleThetaPhiPlot3D lost3d;
        TH1D* angle_regions;
        TH1D* nlost;
        TH1D* energy_reco;
        ParticleThetaPhiPlot matched_pos_after_geo;
        double emin;


        AcceptanceAnalysis(HistogramFactory& factory, const utils::A2SimpleGeometry& geo_, const std::string& name_);
        void Fill(const TParticleList& mctrue, const TParticleList& reconstructed);
        void ShowResult();
    };

    utils::A2SimpleGeometry geo;

    std::list<AcceptanceAnalysis> analyses;

public:
    GeoAcceptance(const std::string& name, OptionsPtr opts);
    virtual ~GeoAcceptance();
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};
}
}
}
