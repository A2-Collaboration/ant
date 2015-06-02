#ifndef DATAOVERVIEW_H
#define DATAOVERVIEW_H

#include "AntPhysics.h"
#include "plot/SmartHist.h"

#include <string>

namespace ant {
namespace analysis {

class DataOverview : public ant::Physics {
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

        void Fill(const ant::Event::Data& dataset);
    };

    OverviewSet reconstructed;
    OverviewSet mctrue;

public:
    DataOverview(const std::string& name="OverviewSet");

    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}

#endif
