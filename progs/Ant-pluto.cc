/**
  *@file Ant-pluto.cc
  *@brief Pluto generator forntend.
  *
  * Two modes available: Single reaction mode and Random Gun
  *
  * Random gun:
  *  Run with --gun to shoot particles of a selected type uniformly
  *  in all directions.
  *
  * Reaction Mode:
  *  Generate physics using pluto by specifying a reaction string.
  *  Run with --pluto --reaction "...."
  *  Where "...." is a Pluto reaction string, for example "p omega [ pi0 g ]":
  *  for omega production followed by a decay into pi0 gamma.
  *
  *  To further decay the particles use --enableBulk. Then all unstable particles
  *  decay according to the branching ratios in the database.
  **/

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"
#include "base/vec/vec3.h"
#include "mc/pluto/PlutoExtensions.h"
#include "mc/pluto/utils/PlutoTID.h"


// pluto++
#include "PParticle.h"
#include "PBeamSmearing.h"
#include "PReaction.h"
#include "PPlutoBulkDecay.h"
#include "PDataBase.h"

// ROOT
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

// Detail
#include "detail/McAction.h"

#include <string>
#include <memory>

using namespace std;
using namespace ant;
using namespace ant::std_ext;

const constexpr auto GeV = 1000.0;
const constexpr auto me  = 0.000510999;  // mass electron [GeV]
/**
 * @brief thetaCrit
 * @return
 */
inline constexpr Double_t thetaCrit(const double Egmax) noexcept {
    return me / Egmax;
}


struct PlutoAction : McAction {
    string   beam;
    string   target;
    string   reaction;
    bool     saveIntermediate;
    bool     enableBulk;
    bool     flatEbeam;
    bool     beamSmearing;
    int      verbosity_level;
    virtual void Run() const override;
};



int main( int argc, char** argv ) {
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-pluto - Pluto single reaction generator for A2 Physics", ' ', "0.1");

    // common options
    auto cmd_numEvents = cmd.add<TCLAP::ValueArg<unsigned>>  ("n", "numEvents", "Number of generated events", true, 0, "unsigned int");
    auto cmd_outfile   = cmd.add<TCLAP::ValueArg<string>>    ("o", "outfile", "Output file", true, "pluto.root", "string");
    auto cmd_Emin      = cmd.add<TCLAP::ValueArg<double>>    ("",  "Emin", "Minimum incident energy [MeV]", false, 0.1, "double [MeV]");
    auto cmd_Emax      = cmd.add<TCLAP::ValueArg<double>>    ("",  "Emax", "Maximum incident energy [MeV]", false, 1.6*GeV, "double [MeV]");
    auto cmd_noTID     = cmd.add<TCLAP::SwitchArg>           ("",  "noTID", "Don't add TID tree for the events", false);
    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>       ("v", "verbose","Verbosity level (0..9)", false, 0,"int");

    // reaction simulation options
    auto cmd_reaction = cmd.add<TCLAP::ValueArg<string>> ("", "reaction", "Pseudo Beam - decay string (reaction string), e.g. 'p pi0 [g g]' for pion photoproduction", true, "", "g p decay string");
    auto cmd_noUnstable = cmd.add<TCLAP::SwitchArg>  ("", "no-unstable", "Don't save unstable/intermediate particles", false);
    auto cmd_noBulk = cmd.add<TCLAP::SwitchArg>      ("", "no-bulk", "Disable bulk decay of particles", false);
    auto cmd_flatEbeam = cmd.add<TCLAP::SwitchArg>   ("", "flatEbeam", "Make tagger spectrum flat (no 1/Ebeam weighting)", false);
    auto cmd_noBeamSmear = cmd.add<TCLAP::SwitchArg> ("", "no-beam-smearing", "Disable angular beam smearing according to Bremsstrahlung", false);


    auto cmd_beam   = cmd.add<TCLAP::ValueArg<string>> ("b", "beam", "Pluto-name for beam-particle type", false, "g", "Pluto-particle name");
    auto cmd_target = cmd.add<TCLAP::ValueArg<string>> ("t", "target", "Pluto-name for target-particle type", false, "p", "Pluto-particle name");


    cmd.parse(argc, argv);

    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    PlutoAction action;
    action.beam             = cmd_beam->getValue();
    action.target           = cmd_target->getValue();
    action.reaction         = cmd_reaction->getValue();
    action.saveIntermediate = !(cmd_noUnstable->getValue());
    action.enableBulk       = !(cmd_noBulk->getValue());
    action.flatEbeam        = cmd_flatEbeam->getValue();
    action.beamSmearing     = !(cmd_noBeamSmear->getValue());
    action.verbosity_level  = cmd_verbose->getValue();


    action.nEvents = cmd_numEvents->getValue();
    action.outfile = cmd_outfile->getValue();
    action.Emin    = cmd_Emin->getValue();
    action.Emax    = cmd_Emax->getValue();

    VLOG(2) << "gRandom is a " << gRandom->ClassName();
    gRandom->SetSeed();  // Initialize ROOT's internal rng. Used for TF1s.
    VLOG(2) << "gRandom initialized";

    action.Run();

    LOG(INFO) << "Simulation finished.";

    // add TID tree for the generated events
    if(!cmd_noTID->isSet()) {
        LOG(INFO) << "Add TID tree to the output file";
        string outfile = action.outfile;
        if(!string_ends_with(outfile, ".root"))
            outfile += ".root";

        mc::pluto::utils::PlutoTID::AddTID(outfile);
    }

    return EXIT_SUCCESS;
}


void PlutoAction::Run() const
{

    VLOG(1) << "Running PlutoGun";
    VLOG(1) << "Number of Events:    " << nEvents;
    VLOG(1) << "Reaction: " << reaction;
    VLOG(1) << "Photon Beam E min: " << Emin << " MeV";
    VLOG(1) << "Photon Beam E max: " << Emax << " MeV";
    VLOG(1) << "Save Intermediate: " << saveIntermediate;
    VLOG(1) << "Enable bulk decays: " << enableBulk;
    VLOG(1) << "Flat Ebeam spectrum (no 1/E weighting): " << flatEbeam;
    VLOG(1) << "Angular smearing of the beam: " << beamSmearing;


#ifdef PLUTO_GLOBAL_H
    pluto_global::verbosity = verbosity_level;
#endif

    PBeamSmearing *smear = new PBeamSmearing(strdup("beam_smear"), strdup("Beam smearing"));
    smear->SetReaction(strdup("g + p"));

    TF1* tagger_spectrum = new TF1("bremsstrahlung", flatEbeam ? "(1.0)" : "(1.0/x)", Emin/GeV, Emax/GeV);
    smear->SetMomentumFunction( tagger_spectrum );

    if (beamSmearing) {
        TF1* theta_smear = new TF1( "angle", "x / ( x*x + [0] )^2", 0.0, 5.0 * thetaCrit(Emax/GeV) );
        theta_smear->SetParameter( 0, thetaCrit(Emax/GeV) * thetaCrit(Emax/GeV) );

        smear->SetAngularSmearing( theta_smear );
    }

    makeDistributionManager()->Add(smear);

    ant::mc::pluto::UpdatePlutoDataBase();

    // remove file ending because pluto attaches a ".root"...
    string outfile_clean(outfile);
    if(string_ends_with(outfile_clean, ".root")) {
        outfile_clean = outfile_clean.substr(0,outfile_clean.size()-5);
    }

    /// @note: not using unique_ptr here because deleting the PReaction crashes. Leeking on purpose here.
    PReaction* Pluto_reaction = new PReaction(Emax/GeV, strdup(beam.c_str()), strdup(target.c_str()),
                                              strdup(reaction.c_str()), strdup(outfile_clean.c_str()),
                                              saveIntermediate, 0, 0, 0);

    if( enableBulk ) {
        PPlutoBulkDecay* p1 = new PPlutoBulkDecay();
        p1->SetRecursiveMode(1);
        p1->SetTauMax(0.001);
        Pluto_reaction->AddBulk(p1);
    }

    Pluto_reaction->Print();

    Pluto_reaction->Loop(nEvents, 0, verbosity_level);

}
