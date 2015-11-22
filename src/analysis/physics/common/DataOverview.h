#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/SmartHist.h"

#include <string>

namespace ant {
namespace analysis {
namespace physics {

class DataOverview : public Physics {
protected:
    class OverviewSet {
    public:
        SmartHist1<int> TaggerChannel;
        SmartHist1<double> PhotonEnergy;
        SmartHist1<double> TaggedTime;

        SmartHist1<int> nParticles;

        SmartHist1<double> CBEnergySum;

        SmartHist1<std::string> ParticleTypes;

        OverviewSet(SmartHistFactory& factory, const std::string& title);

        void Fill(const data::Event::Data& dataset);
    };

    OverviewSet reconstructed;
    OverviewSet mctrue;

public:
    DataOverview(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event);
    void Finish();
    void ShowResult();
};

class TaggerOverview : public Physics {
protected:
    TH1D* nHitsEvent = nullptr;
    TH1D* Channels   = nullptr;
    TH1D* Energies   = nullptr;
    TH1D* Times      = nullptr;

    TH2D* channel_correlation = nullptr;

    enum class Mode { MCTrue, Reconstructed };

    Mode mode = Mode::Reconstructed;

    std::string GetMode() const;

public:
    TaggerOverview(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

}
}
}
