#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/Histogram.h"
#include "analysis/plot/SmartHist.h"
#include "analysis/plot/HistogramFactories.h"
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

        ParticleThetaPhiPlot(SmartHistFactory& factory,
                             const std::string& title,
                             const std::string& name="",
                             const BinSettings& thetabins = BinSettings(180),
                             const BinSettings& phibins = BinSettings(360,-180.0,180.0));
        void Fill(const data::ParticlePtr& p);


        TObject *GetObject();
        void Draw(const std::string &option) const;
    };

    class ParticleThetaPhiPlot3D: public ant::root_drawable_traits {
    public:
        TH3* hist;
        unsigned int n;

        ParticleThetaPhiPlot3D(SmartHistFactory& factory, const std::string& title, const std::string& name="");
        void Fill(const data::ParticlePtr& p);


        TObject *GetObject();
        void Draw(const std::string &option) const;
    };

    class AcceptanceAnalysis {
    public:
        std::string name;
        SmartHistFactory HistFac;
        const utils::A2SimpleGeometry& geo;
        ParticleThetaPhiPlot mctrue_pos;
        ParticleThetaPhiPlot matched_pos;
        ParticleThetaPhiPlot multimatched_pos;
        ParticleThetaPhiPlot lost_pos;
        ParticleThetaPhiPlot lost_pos_zoom;
        ParticleThetaPhiPlot3D lost3d;
        TH1D* angle_regions;
        TH1D* nlost;
        SmartHist1<double> energy_reco;
        ParticleThetaPhiPlot matched_pos_after_geo;
        double emin;


        AcceptanceAnalysis(SmartHistFactory& factory, const utils::A2SimpleGeometry& geo_, const std::string& name_);
        void Fill(const data::ParticleList& mctrue, const data::ParticleList& reconstructed);
        void ShowResult();
    };

    utils::A2SimpleGeometry geo;

    std::list<AcceptanceAnalysis> analyses;

public:
    GeoAcceptance(const std::string& name, PhysOptPtr opts);
    virtual ~GeoAcceptance();
    void ProcessEvent(const data::Event &event);
    void Finish();
    void ShowResult();
};
}
}
}
