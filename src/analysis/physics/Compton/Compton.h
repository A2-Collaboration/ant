// ----------------------- Important -----------------------
// In Compton folder is a README file. Please read it if you
// want to use this code.


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

// Creating a new class called Compton that inherits
// the members of the Physics class. The Physics class is
// defined in "physics.h"
class Compton : public Physics {
public:

// -------------- The Generic Ant Fucntions --------------

    // Constructor declared
    Compton(const std::string& name, OptionsPtr opts);

    // ProcessEvent is a funtion that is used by every physics
    // class. "override" tells the compiler to used this
    // ProcessEvent as default when it runs this class
    // default. TEvent and manager_t are data types that are
    // defined in ant (like int or char or double).
    // "event" is a data structure with all the info for an
    // event. ProcessEvent loops over all the events
    virtual void ProcessEvent(const TEvent& event,
                              manager_t& manager) override;

    // For outputting stuff (like histograms)
    virtual void ShowResult() override;

// ------- Methods specific to the Compton class -------

    // (explinations in cc file)
    bool IsParticleCharged(double veto_energy);

    int IsChargedUncharged(const TCandidateList& candidates);

    double GetMissingMass(const TCandidate& candidate,
                          const LorentzVec target,
                          const LorentzVec incoming);

    double GetCloserMM(const TCandidateList& candidates,
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
//    TH1D* h_WeightedTaggerTime;

    // Preliminary cuts
    // Note: MM means missing mass
//    TH1D* h_MM;
//    TH1D* h_MM1;
//    TH1D* h_MM11;

    // 1 Particle cuts
//    TH1D* h_MM101;
    TH1D* h_MM111;

    // Preliminary 2 particle cuts
//    TH1D* h_MM102;
//    TH1D* h_MM112;
//    TH1D* h_MM1021;

    // Coplanar cuts
//    TH1D* h_MM10201;
//    TH1D* h_MM11201;
//    TH1D* h_MM10211;

    // Opening angle cuts
//    TH1D* h_MM102001;
//    TH1D* h_MM112001;
//    TH1D* h_MM102011;
    TH1D* h_MM112011;
    // Uncharged/Charged cut done before open ang cut
//    TH1D* h_MM112001_switch;
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

    // Stuff for PR cut
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;


// ----- Default values for options at the command line -----

    // Incoming photon energy range cut
    double tagger_energy_low = 0;     // in MeV
    double tagger_energy_high = 2000;

    // Prompt random windows
    std::string PR_windows = "-200,-6,-5,5,6,200";   // in ns

// ----------------- Scalar Counter Objects -----------------

    const std::shared_ptr<TaggerDetector_t> tagger;
    unsigned nchannels = 0;

// ------------------- Other Objects used -------------------

    const double proton_mass = ParticleTypeDatabase::Proton.Mass();

    // Momentum 4 vectors for target (i.e. stationary proton)
    // and incoming photon
    const LorentzVec target_vec = LorentzVec({0.0,0.0,0.0},
                     proton_mass);
    LorentzVec incoming_vec;

    double missing_mass;
    double closer_missing_mass;

};

}}}
