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
  *  Generate physics using pluto by specifiying a reaction string.
  *  Run with --pluto --reaction "...."
  *  Where "...." is a Pluto reaction string, for example "p omega [ pi0 g ]":
  *  for omega production followed by a decay into pi0 gamma.
  *
  *  To further decay the particles use --enableBulk. Then all instable particles
  *  decay according to the database.
  **/

#include "base/CmdLine.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"
#include "base/vec/vec3.h"
#include "simulation/mc/PlutoExtensions.h"


// pluto++
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wvla-extension"
#endif

#include "PParticle.h"
#include "PBeamSmearing.h"
#include "PReaction.h"
#include "PPlutoBulkDecay.h"
#include "PDataBase.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#pragma GCC diagnostic pop

// ROOT
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

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

struct Action {
    unsigned nEvents;
    string   outfile;
    double   Emin;
    double   Emax;
    virtual void Run() const =0;
    virtual ~Action() = default;
};

struct PlutoAction : Action {
    string   reaction;
    bool     saveIntermediate;
    bool     enableBulk;
    virtual void Run() const override;
};

struct GunAction : Action {
    std::vector<string>   particles;
    unsigned nPerEvents;
    double   thetaMin;
    double   thetaMax;
    virtual void Run() const override;
    const static std::map<string,int> available_particles;
};

// these are the allowed particles for the gun with their pluto IDs
const std::map<string,int> GunAction::available_particles = {
    {"p", 14},
    {"g",  1},
    {"e-", 3},
    {"e+", 4},
    {"pi+",8},
    {"pi-",9}
};

/**
 *@breif get a randomly choosen element from a vector
 *@param v the vetor to draw from
 *@return const reference to the element
 */
template <typename T>
const T& getRandomFrom(const std::vector<T>& v) {
    const size_t index = floor(gRandom->Uniform(v.size()));
    return v.at(index);
}

/**
 *@brief get all keys of a map in a vector
 *@param m the map to extract keys from
 *@return vector of copies of the keys
 */
template <typename T, typename U>
std::vector<T> getKeys(const std::map<T,U>& m) {
    vector<T> keys;
    keys.reserve(m.size());
    for(const auto& e : m) {
        keys.emplace_back(e.first);
    }
    return keys;
}

int main( int argc, char** argv ) {
    SetupLogger();

    gRandom->SetSeed(); // Initialize ROOT's internal rng. Used for TF1s.

    TCLAP::CmdLine cmd("Ant-pluto - Pluto single reaction generator for A2 Physics", ' ', "0.1");

    // common options
    auto cmd_numEvents = cmd.add<TCLAP::ValueArg<unsigned>>  ("n", "numEvents", "Number of generated events", true, 0, "unsigned int");
    auto cmd_outfile   = cmd.add<TCLAP::ValueArg<string>>    ("o", "outfile", "Output file", true, "pluto.root", "string");
    auto cmd_Emin      = cmd.add<TCLAP::ValueArg<double>>    ("",  "Emin", "Minimal incident energy [MeV]", false, 0.0, "double [MeV]");
    auto cmd_Emax      = cmd.add<TCLAP::ValueArg<double>>    ("",  "Emax", "Maximal incident energy [MeV]", false, 1.6*GeV, "double [MeV]");
    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>       ("v", "verbose","Verbosity level (0..9)", false, 0,"int");

    // reaction simulation options
//    auto cmd_usePluto = cmd.add<TCLAP::SwitchArg>        ("", "pluto", "Use pluto", false);
    auto cmd_reaction = cmd.add<TCLAP::ValueArg<string>> ("", "reaction", "Reaction string in PLUTO notation", false, "", "string (Pluto Rection)");
    auto cmd_noUnstable = cmd.add<TCLAP::SwitchArg>("", "no-unstable", "Don't unstable particles", false);
    auto cmd_noBulk = cmd.add<TCLAP::SwitchArg>      ("", "no-bulk", "Disable bulk decay of particles", false);

    // random gun options
    auto cmd_useGun = cmd.add<TCLAP::SwitchArg>        ("", "gun", "Use random gun", false);

    TCLAP::ValuesConstraintExtra<vector<string>> allowed_particles(getKeys(GunAction::available_particles));
    auto cmd_randomparticles = cmd.add<TCLAP::MultiArg<string>>  ("p", "particle", "Particle type to shoot", false, &allowed_particles);
    auto cmd_numParticles   = cmd.add<TCLAP::ValueArg<unsigned>> ("N", "particles-event", "Paricles per event", false, 1, "unsigned int");
    auto cmd_thetaMin       = cmd.add<TCLAP::ValueArg<double>>   ("",  "theta-min", "Minimal theta angle [deg]", false,   0.0, "double [deg]");
    auto cmd_thetaMax       = cmd.add<TCLAP::ValueArg<double>>   ("",  "theta-max", "Maximal theta angle [deg]", false, 180.0, "double [deg]");


    cmd.parse(argc, argv);

    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    if(cmd_useGun->isSet() && cmd_reaction->isSet()) {
        LOG(ERROR) << "Can't use Pluto and random gun at the same time!";
        return EXIT_FAILURE;
    }

    if(!cmd_useGun->isSet() && !cmd_reaction->isSet()) {
        LOG(ERROR) << "Require one of --pluto , --gun";
        return EXIT_FAILURE;
    }

    unique_ptr<Action> action;

    if(cmd_useGun->isSet()) {
        unique_ptr<GunAction> gunaction = std_ext::make_unique<GunAction>();
        if(cmd_randomparticles->getValue().empty()) {
            LOG(ERROR) << "Need particles for random gun! (--particle xxx)";
            return EXIT_FAILURE;
        }
        gunaction->particles  = cmd_randomparticles->getValue();
        gunaction->nPerEvents = cmd_numParticles->getValue();
        gunaction->thetaMin   = degree_to_radian(cmd_thetaMin->getValue());
        gunaction->thetaMax   = degree_to_radian(cmd_thetaMax->getValue());
        action = move(gunaction);
    } else if(cmd_reaction->isSet()){
        auto plutoaction              = std_ext::make_unique<PlutoAction>();
        plutoaction->reaction         = cmd_reaction->getValue();
        plutoaction->saveIntermediate = !(cmd_noUnstable->getValue());
        plutoaction->enableBulk       = !(cmd_noBulk->getValue());
        action = move(plutoaction);
    }

    assert(action);

    action->nEvents = cmd_numEvents->getValue();
    action->outfile = cmd_outfile->getValue();
    action->Emin    = cmd_Emin->getValue();
    action->Emax    = cmd_Emax->getValue();

    VLOG(2) << "gRandom is a " << gRandom->ClassName();
    gRandom->SetSeed();
    VLOG(2) << "gRandom initialized";

    action->Run();

    cout << "Simulation finished." << endl;

    // Do not delete the reaction, otherwise: infinite loop somewhere in ROOT...
    //delete reactrion;

    return 0;
}


void PlutoAction::Run() const
{

    VLOG(1) << "Running PlutoGun";
    VLOG(1) << "number of Events:    " << nEvents;
    VLOG(1) << "Reaction: " << reaction;
    VLOG(1) << "Photon Beam E min: " << Emin << " MeV";
    VLOG(1) << "Photon Beam E max: " << Emax << " MeV";
    VLOG(1) << "Save Intermediate: " << saveIntermediate;
    VLOG(1) << "Enable bulk decays: " << enableBulk;

    PBeamSmearing *smear = new PBeamSmearing(strdup("beam_smear"), strdup("Beam smearing"));
    smear->SetReaction(strdup("g + p"));

    TF1* tagger_spectrum = new TF1("bremsstrahlung","(1.0/x)", Emin/GeV, Emax/GeV);
    smear->SetMomentumFunction( tagger_spectrum );

    TF1* theta_smear = new TF1( "angle", "x / ( x*x + [0] )^2", 0.0, 5.0 * thetaCrit(Emax/GeV) );
    theta_smear->SetParameter( 0, thetaCrit(Emax/GeV) * thetaCrit(Emax/GeV) );

    smear->SetAngularSmearing( theta_smear );

    makeDistributionManager()->Add(smear);

    ant::simulation::mc::UpdatePluteDataBase();

    // remove file ending because pluto attaches a ".root"...
    string outfile_clean(outfile);
    if(string_ends_with(outfile_clean, ".root")) {
        outfile_clean = outfile_clean.substr(0,outfile_clean.size()-5);
    }

    /// @note: not using unique_ptr here because deleting the PReaction crashes. Leeking on purpose here.
    PReaction* Pluto_reaction = new PReaction(Emax/GeV, strdup("g"), strdup("p"),
                                        strdup(reaction.c_str()), strdup(outfile_clean.c_str()),
                                        saveIntermediate, 0, 0, 0);

    if( enableBulk ) {
        PPlutoBulkDecay* p1 = new PPlutoBulkDecay();
        p1->SetRecursiveMode(1);
        p1->SetTauMax(0.001);
        Pluto_reaction->AddBulk(p1);
    }

    Pluto_reaction->Print();   //The "Print()" statement is optional

    Pluto_reaction->Loop(nEvents);

    //Pluto_reaction->Close(); // no using -> crashes

}

void GunAction::Run() const
{
    VLOG(1) << "Running RandomGun";
    VLOG(1) << "number of Events:    " << nEvents;
    VLOG(1) << "number of particles: " << nPerEvents;
    VLOG(1) << "Particles: " << particles;
    VLOG(1) << "E min: " << Emin << " MeV";
    VLOG(1) << "E max: " << Emax << " MeV";
    VLOG(1) << "Theta min: " << radian_to_degree(thetaMin) << " degree";
    VLOG(1) << "Theta max: " << radian_to_degree(thetaMax) << " degree";

    auto db = makeStaticData();

    WrapTFileOutput file(outfile, WrapTFileOutput::mode_t::recreate, false);

    TTree* tree = file.CreateInside<TTree>("data","Random Particles");

    TClonesArray* particles_array = new TClonesArray("PParticle", nPerEvents);
    particles_array->SetOwner(kTRUE);

    tree->Branch("Particles", particles_array);

    std::vector<int> ids;
    ids.reserve(particles.size());
    for(const auto& p : particles) {
        ids.push_back(available_particles.at(p));
    }

    for( unsigned i=0; i< nEvents; ++i ) {

        particles_array->Clear();

        for( unsigned j=0; j<nPerEvents; ++j ) {

            const int pID = ids.size()==1 ? ids[0] : getRandomFrom(ids);

            const double m = db->GetParticleMass(pID);
            const double E = gRandom->Uniform(Emax/GeV)+m;

            const double p = sqrt(E*E - m*m);

            ant::vec3 dir;
            do {
                gRandom->Sphere(dir.x, dir.y, dir.z, p);
            } while (dir.Theta() > thetaMax || dir.Theta() < thetaMin);

            PParticle* part = new PParticle( pID, dir);

            (*particles_array)[j] = part;
        }

        tree->Fill();

    }
}
