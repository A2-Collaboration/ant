#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/SmartHist.h"

#include <string>

namespace ant {

class TH2CB;

namespace analysis {
namespace physics {


/**
 * @brief Base class for DataOverview physics
 *
 * Options:
 *    Mode = [ MCTrue | Reconstructed ] : Analyse mc true or reconstructed data branch.
 *         The default is "Reconstructed".
 */
class DataOverviewBase : public Physics {
protected:

    enum class Mode { MCTrue, Reconstructed };

    Mode mode = Mode::Reconstructed;

    /**
     * @brief Get a string representation of the current mode
     * @return "MCTrue" or "Reconstructed"
     */
    std::string GetMode() const;

    const data::Event::Data& GetBranch(const data::Event& event) const;

public:
    DataOverviewBase(const std::string& name, PhysOptPtr opts);
    virtual ~DataOverviewBase();
};

/**
 * @brief Physics class to show Tagger related infos
 */
class TaggerOverview : public DataOverviewBase {
protected:
    TH1D* nHitsEvent = nullptr;
    TH1D* Channels   = nullptr;
    TH1D* Energies   = nullptr;
    TH1D* Times      = nullptr;

    TH2D* channel_correlation = nullptr;

public:
    TaggerOverview(const std::string& name, PhysOptPtr opts);
    virtual ~TaggerOverview();

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

/**
  * @brief Physics class to show Trigger and general event infos
  */
class TriggerOverview : public DataOverviewBase {
protected:
    TH1D* CBESum       = nullptr;
    TH1D* Multiplicity = nullptr;
    TH1D* nErrorsEvent = nullptr;
    ///@todo Add histograms for Error Codes and ModuleIDs ?

    TH2CB* oneclusterevents = nullptr;

public:
    TriggerOverview(const std::string& name, PhysOptPtr opts);
    virtual ~TriggerOverview();

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

/**
 * @brief The TargetOverview class shows target information
 */
class TargetOverview : public DataOverviewBase {
protected:
    TH2D* VertexXY;
    TH1D* VertexZ;

public:
    TargetOverview(const std::string& name, PhysOptPtr opts);
    virtual ~TargetOverview();

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

/**
 * @brief The ParticleOverview class shows particle type occurences
 */
class ParticleOverview: public DataOverviewBase {
protected:
    TH1D* nParticles    = nullptr;
    TH1D* particleTypes = nullptr;

    std::map<const ParticleTypeDatabase::Type*, TH1D*> nType;

    static void SetBinLabels(TH1D* hist, const ParticleTypeDatabase::TypeList_t& types);

public:
    ParticleOverview(const std::string& name, PhysOptPtr opts);
    virtual ~ParticleOverview();

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};



}
}
}
