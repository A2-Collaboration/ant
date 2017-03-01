#include "base/CmdLine.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"
#include "base/vec/vec3.h"
#include "base/ParticleType.h"

// pluto
#include "PParticle.h"
#include "PDataBase.h"

// ROOT
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

// detail
#include "detail/McAction.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;

const constexpr auto GeV = 1000.0;

using particle_type_list_t = std::vector<const ParticleTypeDatabase::Type*>;

// these are the allowed particles for the gun with their pluto IDs
const particle_type_list_t available_particles = [] () {
    particle_type_list_t v;

    for(auto& type : ParticleTypeDatabase()) {
        v.push_back(addressof(type));
    }

    return v;
}();

std::vector<std::string> get_available_particle_names() {
    std::vector<std::string> names;
    for(auto p : available_particles)
        names.push_back(p->PlutoName());
    return names;
}

struct GunAction : McAction {
    std::vector<string>   particles;
    double   thetaMin;
    double   thetaMax;
    virtual void Run() const override;
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



int main( int argc, char** argv ) {
    SetupLogger();

    gRandom->SetSeed(); // Initialize ROOT's internal rng. Used for TF1s.

    TCLAP::CmdLine cmd("Ant-mcgun - Simple particle gun.", ' ', "0.1");

    // random gun options
    TCLAP::ValuesConstraintExtra<vector<string>> allowed_particles(get_available_particle_names());
    auto cmd_randomparticles = cmd.add<TCLAP::MultiArg<string>>  ("p", "particle", "Particle type to shoot", true, &allowed_particles);

    auto cmd_thetaMin       = cmd.add<TCLAP::ValueArg<double>>   ("",  "theta-min", "Minimal theta angle [deg] for first particle in an event.", false,   0.0, "double [deg]");
    auto cmd_thetaMax       = cmd.add<TCLAP::ValueArg<double>>   ("",  "theta-max", "Maximal theta angle [deg] for first particle in an event.", false, 180.0, "double [deg]");

    // common options
    auto cmd_numEvents = cmd.add<TCLAP::ValueArg<unsigned>>  ("n", "numEvents", "Number of generated events", true, 0, "unsigned int");
    auto cmd_outfile   = cmd.add<TCLAP::ValueArg<string>>    ("o", "outfile", "Output file", true, "pluto.root", "string");
    auto cmd_Emin      = cmd.add<TCLAP::ValueArg<double>>    ("",  "Emin", "Minimal incident energy [MeV]", false, 0.0, "double [MeV]");
    auto cmd_Emax      = cmd.add<TCLAP::ValueArg<double>>    ("",  "Emax", "Maximal incident energy [MeV]", false, 1.6*GeV, "double [MeV]");
    auto cmd_OpenAngle = cmd.add<TCLAP::ValueArg<double>>    ("",  "OpeningAngle", "Maximal opening angle to first produced particle in an event.", false, 180.0, "double [deg]");

    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>       ("v", "verbose","Verbosity level (0..9)", false, 0,"int");




    cmd.parse(argc, argv);

    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }



    GunAction action;

    action.particles  = cmd_randomparticles->getValue();
    action.thetaMin   = degree_to_radian(cmd_thetaMin->getValue());
    action.thetaMax   = degree_to_radian(cmd_thetaMax->getValue());


    action.nEvents = cmd_numEvents->getValue();
    action.outfile = cmd_outfile->getValue();
    action.Emin    = cmd_Emin->getValue();
    action.Emax    = cmd_Emax->getValue();

    VLOG(2) << "gRandom is a " << gRandom->ClassName();
    gRandom->SetSeed();
    VLOG(2) << "gRandom initialized";

    action.Run();

    cout << "Simulation finished." << endl;

    return 0;

}

void GunAction::Run() const
{
    auto nParticles = particles.size();

    VLOG(1) << "Running RandomGun";
    VLOG(1) << "number of Events:    " << nEvents;
    VLOG(1) << "number of particles: " << nParticles;
    VLOG(1) << "Particles: " << particles;
    VLOG(1) << "E min: " << Emin << " MeV";
    VLOG(1) << "E max: " << Emax << " MeV";
    VLOG(1) << "Theta min: " << radian_to_degree(thetaMin) << " degree";
    VLOG(1) << "Theta max: " << radian_to_degree(thetaMax) << " degree";


    WrapTFileOutput file(outfile, WrapTFileOutput::mode_t::recreate, false);

    TTree* tree = file.CreateInside<TTree>("data","Random Particles");

    TClonesArray* particles_array = new TClonesArray("PParticle", int(particles.size()));
    particles_array->SetOwner(kTRUE);

    tree->Branch("Particles", particles_array);

    std::vector<int> ids;
    ids.reserve(particles.size());
    for(const auto& p : particles) {
//        ids.push_back(available_particles.at(p));
    }

    for( unsigned evt=0; evt< nEvents; ++evt ) {

        particles_array->Clear();

        for( unsigned iParticle=0; iParticle<nParticles; ++iParticle ) {

            const int pID = ids.size()==1 ? ids[0] : getRandomFrom(ids);

            const double m = 0;

            const double E = gRandom->Uniform(Emax/GeV)+m;

            const double p = sqrt(E*E - m*m);

            ant::vec3 dir;
            do {
                gRandom->Sphere(dir.x, dir.y, dir.z, p);
            } while (dir.Theta() > thetaMax || dir.Theta() < thetaMin);

            PParticle* part = new PParticle( pID, dir);

            (*particles_array)[iParticle] = part;
        }

        tree->Fill();

    }
}
