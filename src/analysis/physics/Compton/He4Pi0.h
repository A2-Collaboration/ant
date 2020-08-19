// If this file is #include more than once, this command tells
// the compiler to read this file only once
#pragma once

#include "physics/Physics.h"

// To subtact out random tagger hits
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

// To get stuff at the command line
#include "base/Logger.h"

// Scaler Counts
#include "slowcontrol/SlowControlVariables.h"

// A heirarchy of namespaces which generally resembles the
// folder heirarchy
namespace ant {
namespace analysis {
namespace physics {

// Creating a new class called He4Pi0 that inherits
// the members of the Physics class. The Physics class is
// defined in "physics.h"
class He4Pi0 : public Physics {
public:

// -------------- The Generic Ant Fucntions --------------

    // Constructor declared
    He4Pi0(const std::string& name, OptionsPtr opts);

    // ProcessEvent is a funtion that is used by every physics
    // class. "override" tells the compiler to used this
    // ProcessEvent as default when it runs this class
    // default. TEvent and manager_t are data types that are
    // defined in ant (like int or char or double).
    // "event" is a data structure with all the info for an
    // event. ProcessEvent loops over all the events
    virtual void ProcessEvent(const TEvent& event,
                              manager_t& manager) override;

    // For last minute modifications and such that you don't
    // need to do when looping through every event
    virtual void Finish() override;

    // For outputting stuff (like histograms)
    virtual void ShowResult() override;

// -------- Methods specific to the Compton class --------

    // (explinations in cc file)
    bool IsParticleCharged(double veto_energy);

    bool IsTwoPhotons(const TCandidateList& candidates);

    int IsChargedUncharged(const TCandidateList& candidates);

    LorentzVec GetLorentzVec(const TCandidate& candidate);

    double GetPi0MissingMass(const TCandidate& front_photon,
                             const TCandidate& back_photon);

    LorentzVec GetPi0Vec(const TCandidate& front_photon,
                         const TCandidate& back_photon);

    double GetHe4MissingEnergy(const TCandidate& front_photon,
                               const TCandidate& back_photon,
                               const LorentzVec target,
                               const LorentzVec incoming);

    double GetPi0MissingEnergy(const LorentzVec pi0,
                               const LorentzVec target,
                               const LorentzVec photon);

    double GetMissingMass(const TCandidate& candidate,
                          const LorentzVec target,
                          const LorentzVec incoming);

    bool IsCoplanar(const TCandidateList& candidates);

    int IsOpeningAngle(const TCandidateList& candidates,
                       const LorentzVec target,
                       const LorentzVec incoming);

    bool IsOpeningAngle2(const TCandidateList& candidates,
                         const LorentzVec target,
                         const LorentzVec incoming,
                         const int IsChargedUncharged_output);

// ------------------- Other Methods -------------------

    void PlotCounts();

private:

// --------- Where all the histograms are declared ---------

    // TH1 is a root command for making histograms. The D at
    // the end stands for double and it indicates what type
    // the counts in each bin will be.

    // Tagger hits with weights applied (to check if PR
    // windows were chosen well)
    TH1D* h_WeightedTaggerTime;

    // 1 Particle cuts
    TH1D* h_MM111;

    // Opening angle cuts
    TH1D* h_MM112011;
    TH1D* h_MM112011_switch;

    // 3D Plots
    TH3D* h3D_MM111;
    TH3D* h3D_MM112011;
    TH3D* h3D_MM112011_switch;

    // 3D Plot Projections
    TH1D* h3D_MM111_projX;
    TH1D* h3D_MM112011_projX;
    TH1D* h3D_MM112011_switch_projX;

    // Scalar Counter
    TH1D* h_ScalarCounts;
    TH1D* h_ScalarCountsEnergy;

    // Pi0 Production Histograms
    TH1D* h_MMpi0;
    TH1D* h_MMhe4;
    TH1D* h_MEpi0;

    // Pi0 Production for 3 Particle Events
    TH1D* h_MMpi0_3;
    TH1D* h_MMhe4_3;
    TH1D* h_MEpi0_3;

    // 3D Pi0 ME Histogram
    TH3D* h3D_MEpi0;
    TH3D* h3D_MEpi0_convert;

    // Stuff for PR cut
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

// ----- Default values for options at the command line -----

    // Incoming photon energy range cut
    double tagger_energy_low = 0;     // in MeV
    double tagger_energy_high = 2000;

    // Prompt random windows
    std::string PR_windows = "-200,-6,-5,5,6,200";   // in ns

    // Target
    std::string target_type = "he4";

// ----------------- Scalar Counter Objects -----------------

    const std::shared_ptr<TaggerDetector_t> tagger;
    unsigned seenScalerBlocks = 0;
    unsigned nchannels = 0;

// ------------------ Useful Known Masses ------------------

    const double proton_mass = ParticleTypeDatabase::Proton.Mass();
    const double He4_mass = 3727.84; // not defined in particle database yet
    const double pi0_mass = ParticleTypeDatabase::Pi0.Mass();
    double target_mass;   // eventually want to just set target mass using setup file

// ------------------- Other Objects used -------------------

    // const LorentzVec target_vec = LorentzVec({0.0,0.0,0.0},proton_mass);
    LorentzVec target_vec;
    LorentzVec incoming_vec;
    LorentzVec pi0_vec;

    double missing_mass;
    double pi0_missing_mass;
    double missing_energy;
    double pi0_E_miss;


};

}}}
