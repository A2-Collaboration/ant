#ifndef GEOACCEPTANCE_H
#define GEOACCEPTANCE_H

#include "AntPhysics.h"
#include "plot/Histogram.h"
#include "plot/SmartHist.h"
#include "plot/HistogramFactories.h"
#include "plot/root_draw.h"

#include "A2GeoAcceptance.h"

class TH3;
class TH1D;

namespace ant {
namespace analysis {



class GeoAcceptance: public ant::Physics {

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
        void Fill(const ParticlePtr& p);


        TObject *GetObject();
        void Draw(const std::string &option) const;
    };

    class ParticleThetaPhiPlot3D: public ant::root_drawable_traits {
    public:
        TH3* hist;
        unsigned int n;

        ParticleThetaPhiPlot3D(SmartHistFactory& factory, const std::string& title, const std::string& name="");
        void Fill(const ParticlePtr& p);


        TObject *GetObject();
        void Draw(const std::string &option) const;
    };

    class AcceptanceAnalysis {
    public:
        std::string name;
        SmartHistFactory HistFac;
        const A2SimpleGeometry& geo;
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


        AcceptanceAnalysis(SmartHistFactory& factory, const A2SimpleGeometry& geo_, const std::string& name_);
        void Fill(const ParticleList& mctrue, const ParticleList& reconstructed);
        void ShowResult();
    };

    A2SimpleGeometry geo;

    std::list<AcceptanceAnalysis> analyses;

public:
    GeoAcceptance(const std::string& name="GeoAcceptance", const mev_t energy_scale=1000.0);
    virtual ~GeoAcceptance();
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};
}
}

#endif

