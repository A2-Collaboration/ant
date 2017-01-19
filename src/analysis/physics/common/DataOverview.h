#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"

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

    const TEventData& GetBranch(const TEvent& event) const;

public:
    DataOverviewBase(const std::string& name, OptionsPtr opts);
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
    TaggerOverview(const std::string& name, OptionsPtr opts);
    virtual ~TaggerOverview();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

/**
  * @brief Physics class to show Trigger and general event infos
  */
class TriggerOverview : public DataOverviewBase {
protected:
    TH1D* CBESum       = nullptr;
    TH1D* Multiplicity = nullptr;
    TH1D* nErrorsEvent = nullptr;
    TH1D* CBTiming     = nullptr;
    ///@todo Add histograms for Error Codes and ModuleIDs ?

    TH2D* CBESum_perCh = nullptr;
    TH2D* E_perCh = nullptr;

    struct Tree : WrapTTree {
        ADD_BRANCH_T(double, CBESum)
        ADD_BRANCH_T(double, TaggE)
        ADD_BRANCH_T(double, TaggT)
        ADD_BRANCH_T(TLorentzVector, true_proton)
    };

    TTree* t = nullptr;
    Tree tree;

public:
    TriggerOverview(const std::string& name, OptionsPtr opts);
    virtual ~TriggerOverview();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

/**
 * @brief The TargetOverview class shows target information
 */
class TargetOverview : public DataOverviewBase {
protected:
    TH2D* VertexXY;
    TH1D* VertexZ;

public:
    TargetOverview(const std::string& name, OptionsPtr opts);
    virtual ~TargetOverview();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
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
    ParticleOverview(const std::string& name, OptionsPtr opts);
    virtual ~ParticleOverview();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};



}
}
}
