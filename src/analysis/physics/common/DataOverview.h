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

}
}
}
