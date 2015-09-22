#include "base/CmdLine.h"

// pluto++
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "PParticle.h"
#include "PBeamSmearing.h"
#include "PReaction.h"
#include "PPlutoBulkDecay.h"
#include "PDataBase.h"
#pragma GCC diagnostic pop

// ROOT
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

// Options (default values)
static Double_t Egmax = 1.604;  // beam energy, max Tagg Energy [GeV]
static Double_t Egmin = 1.450;  // beam energy, min Tagg Energy [GeV]
static bool enableBulk = true;
static bool saveIntermediate = true;
static Int_t numEvents = 1000;
static string reactionString = "p omega";
static string prefix = "pluto_";
static string suffix = "00";
static string outfile = "";

static const Double_t me = 0.000510999;  // mass electron [GeV]

/**
 * @brief thetaCrit
 * @return
 */
inline Double_t thetaCrit() {
    return me / Egmax;
}

std::string generateFilename() {
    stringstream s;
    s << prefix;

    string fixedName(reactionString);
    replace( fixedName.begin(), fixedName.end(), ' ', '_');
    s << "-" << fixedName;

    if(enableBulk)
        s << "-bulk";

    s << "-" << suffix;

    return s.str();
}

void ReadCmdline(int argc, char** argv ) {
    TCLAP::CmdLine cmd("Ant-pluto - Pluto single reaction generator for A2 Physics", ' ', "0.1");

    auto cmd_numEvents = cmd.add<TCLAP::ValueArg<int>>("N", "numEvents", "Number of generated events", false, numEvents, "nEvents");
    auto cmd_saveIntermediate = cmd.add<TCLAP::SwitchArg>("", "saveIntermediate", "Save intermediate particles", saveIntermediate);
    auto cmd_enableBulk = cmd.add<TCLAP::SwitchArg>("", "enableBulk", "Enable bulk decay of particles", enableBulk);
    auto cmd_Emin = cmd.add<TCLAP::ValueArg<double>> ("", "Emin", "Minimal photon energy [GeV]", false, Egmin, "GeV");
    auto cmd_Emax = cmd.add<TCLAP::ValueArg<double>> ("", "Emax", "Maximal photon energy [GeV]", false, Egmax, "GeV");
    auto cmd_reaction = cmd.add<TCLAP::ValueArg<string>> ("", "reaction", "Reaction string in PLUTO notation", true, reactionString, "reaction");
    auto cmd_prefix = cmd.add<TCLAP::ValueArg<string>> ("", "prefix", "Prefix for output file", false, prefix, "prefix");
    auto cmd_suffix = cmd.add<TCLAP::ValueArg<string>> ("", "suffix", "Suffix for output file", false, suffix, "suffix");
    auto cmd_outfile = cmd.add<TCLAP::ValueArg<string>> ("", "outfile", "Output file", false, outfile, "outfile");

    cmd.parse(argc, argv);

    numEvents = cmd_numEvents->getValue();
    saveIntermediate = cmd_saveIntermediate->getValue();
    enableBulk = cmd_enableBulk->getValue();
    Egmax = cmd_Emax->getValue();
    Egmin = cmd_Emin->getValue();
    reactionString = cmd_reaction->getValue();
    prefix = cmd_prefix->getValue();
    suffix = cmd_suffix->getValue();
    outfile = cmd_outfile->getValue();
}

void PrintConfig() {
    cout << "\n\n Simulating " << numEvents << " events:\n\n";
    cout << "  Reaction:  g p --> " << reactionString << "\n\n";
    cout << "  Photon Engery : " << Egmin << " to " << Egmax << " GeV\n\n";
    cout << "  Saving to " << generateFilename() << "\n\n";
    cout << "  saveIntermediate particles: ";
    if( saveIntermediate )
        cout << "yes";
    else
        cout << "no";
    cout << "\n\n";
    cout << "  enable Bulk decay: ";
    if( enableBulk )
        cout << "yes";
    else
        cout << "no";
    cout << "\n\n";

}

int main( int argc, char** argv ) {

    gRandom->SetSeed(); // Initialize ROOT's internal rng. Used for TF1s.

    ReadCmdline( argc, argv );
    PrintConfig();

    PBeamSmearing *smear = new PBeamSmearing(strdup("beam_smear"), strdup("Beam smearing"));
    smear->SetReaction(strdup("g + p"));

    TF1* tagger_spectrum = new TF1("bremsstrahlung","(1.0/x)", Egmin, Egmax);
    smear->SetMomentumFunction( tagger_spectrum );

    TF1* theta_smear = new TF1( "angle", "x / ( x*x + [0] )^2", 0.0, 5.0 * thetaCrit() );
    theta_smear->SetParameter( 0, thetaCrit() * thetaCrit() );

    smear->SetAngularSmearing( theta_smear );

    makeDistributionManager()->Add(smear);

    PStaticData* sdata = makeStaticData();

    // additional omega decays (from PDG Booklet 2014)
    sdata->AddDecay("w --> eta + g",        "w", "eta,g",           4.6E-4);
    sdata->AddDecay("w --> g + g + g",      "w", "g,g,g",           1.9E-4);    // upper limit
    sdata->AddDecay("w --> pi0 e+ e-",      "w", "pi0,e+,e-",       7.7E-4);
    sdata->AddDecay("w --> pi0 mu+ mu-",    "w", "pi0,mu+,mu-",     1.3E-4);
    sdata->AddDecay("w --> pi+ pi- pi0 pi0","w", "pi+,pi-,pi0,pi0", 2.0E-4);    // upper limit
    sdata->AddDecay("w --> pi+ pi- pi+ pi-","w", "pi+,pi-,pi+,pi-", 1.0E-3);    // upper limit
    sdata->AddDecay("w --> pi0 pi0 g",      "w", "pi0,pi0,g",       6.6E-5);
    sdata->AddDecay("w --> eta pi0 g",      "w", "eta,pi0,g",       3.3E-5);    // upper limit

    // Charge conjucation violating modes
    sdata->AddDecay("w --> eta pi0",        "w", "eta,pi0",         2.1E-4);    // upper limit
    sdata->AddDecay("w --> pi0 pi0",        "w", "pi0,pi0",         2.1E-4);    // upper limit
    sdata->AddDecay("w --> pi0 pi0 pi0",    "w", "pi0,pi0,pi0",     2.3E-4);    // upper limit

    if( outfile.empty() )
        outfile = generateFilename();

    cout << "filename " << outfile << endl;

    PReaction* reaction = new PReaction(Egmax, strdup("g"), strdup("p"),
                                        strdup(reactionString.c_str()), strdup(outfile.c_str()),
                                        saveIntermediate, 0, 0, 0);

    if( enableBulk ) {
        PPlutoBulkDecay* p1 = new PPlutoBulkDecay();
        p1->SetRecursiveMode(1);
        p1->SetTauMax(0.001);
        reaction->AddBulk(p1);
    }

    reaction->Print();   //The "Print()" statement is optional

    reaction->Loop(numEvents);

    cout << "Simulation finished." << endl;

    // Do not delete the reaction, otherwise: infinite loop somewhere in ROOT...
    //delete reactrion;

    return 0;
}

