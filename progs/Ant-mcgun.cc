#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"
#include "base/vec/vec3.h"
#include "base/ParticleType.h"
#include "tree/TParticle.h"

#include "mc/pluto/utils/PlutoTID.h"


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

const particle_type_list_t getAntParticles(const vector<string>& ss) {
    particle_type_list_t v;
    for ( auto p: ss)
    {
        for(auto& type : ParticleTypeDatabase()) {
            try{
                if (type.PlutoName() == p)
                {
                    v.push_back(addressof(type));
                    break;
                }
            }
            catch(exception&){}
        }
    }
    return v;
}

std::vector<std::string> get_available_particle_names() {
    std::vector<std::string> names;
    for(auto p : available_particles)
    {
        try{
            auto pName = p->PlutoName();
            if ( pName == "" )
                continue;
            names.push_back(pName);

        }
        catch(exception&) {};
    }
    return names;
}

struct GunAction : McAction {
    particle_type_list_t   particles;
    double   thetaMin;
    double   thetaMax;
    bool     flatTheta;
    double   openAngle;

    vec3 getRandomDir() const
    {
        ant::vec3 dir;
        if (!flatTheta) {
            do {
                gRandom->Sphere(dir.x, dir.y, dir.z,1.0);
            } while (dir.Theta() > thetaMax || dir.Theta() < thetaMin);
        } else {
            const double theta = thetaMin + gRandom->Uniform(thetaMax - thetaMin);
            const double phi = gRandom->Uniform(2*M_PI);
            dir = vec3::RThetaPhi(1., theta, phi);
        }
        return dir;
    }

    double getRandomMomentum(const ParticleTypeDatabase::Type* type) const
    {
        const auto m = type->Mass() / GeV;
        const auto E = (Emin  + gRandom->Uniform( Emax - Emin))/ GeV + m;
        const auto p = sqrt(E*E - m*m);
        return p;
    }

    vec3 getRandomDirInCone(const vec3 center) const
    {
        ant::vec3 dir;
        do {
            gRandom->Sphere(dir.x, dir.y, dir.z,1.0);
        } while (center.Angle(dir) > openAngle);
        return dir;
    }


    virtual void Run() const override;
};



int main( int argc, char** argv ) {
    SetupLogger();

    gRandom->SetSeed(); // Initialize ROOT's internal rng. Used for TF1s.

    TCLAP::CmdLine cmd("Ant-mcgun - Simple particle gun.", ' ', "0.1");

    // random gun options
    TCLAP::ValuesConstraintExtra<vector<string>> allowed_particles(get_available_particle_names());
    auto cmd_particles = cmd.add<TCLAP::MultiArg<string>>      ("p", "particle",     "Particle type to shoot", true, &allowed_particles);

    auto cmd_thetaMin       = cmd.add<TCLAP::ValueArg<double>> ("",  "thetaMin",     "Minimal theta angle [deg] for first particle in an event.", false,   0.0, "double [deg]");
    auto cmd_thetaMax       = cmd.add<TCLAP::ValueArg<double>> ("",  "thetaMax",     "Maximal theta angle [deg] for first particle in an event.", false, 180.0, "double [deg]");
    auto cmd_flatTheta      = cmd.add<TCLAP::SwitchArg>        ("",  "flatTheta",    "Generate a flat theta distribution instead of a random vector on a sphere.", false);

    // common options
    auto cmd_numEvents = cmd.add<TCLAP::ValueArg<unsigned>>    ("n", "numEvents",    "Number of generated events", true, 0, "unsigned int");
    auto cmd_outfile   = cmd.add<TCLAP::ValueArg<string>>      ("o", "outfile",      "Output file", true, "pluto.root", "string");
    auto cmd_Emin      = cmd.add<TCLAP::ValueArg<double>>      ("",  "Emin",         "Minimal incident energy [MeV]", false, 0.0, "double [MeV]");
    auto cmd_Emax      = cmd.add<TCLAP::ValueArg<double>>      ("",  "Emax",         "Maximal incident energy [MeV]", false, 1.6*GeV, "double [MeV]");
    auto cmd_OpenAngle = cmd.add<TCLAP::ValueArg<double>>      ("",  "openingAngle", "Maximal opening angle to first produced particle in an event.", false, std::numeric_limits<double>::quiet_NaN(), "double [deg]");

    auto cmd_noTID     = cmd.add<TCLAP::SwitchArg>             ("",  "noTID",        "Don't add TID tree for the events", false);
    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>         ("v", "verbose",      "Verbosity level (0..9)", false, 0,"int");


    cmd.parse(argc, argv);

    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    GunAction action;

    action.particles  = getAntParticles(cmd_particles->getValue());
    action.thetaMin   = degree_to_radian(cmd_thetaMin->getValue());
    action.thetaMax   = degree_to_radian(cmd_thetaMax->getValue());
    action.openAngle  = degree_to_radian(cmd_OpenAngle->getValue());
    action.flatTheta  = cmd_flatTheta->isSet();


    action.nEvents = cmd_numEvents->getValue();
    action.outfile = cmd_outfile->getValue();
    action.Emin    = cmd_Emin->getValue();
    action.Emax    = cmd_Emax->getValue();

    VLOG(2) << "gRandom is a " << gRandom->ClassName();
    gRandom->SetSeed();
    VLOG(2) << "gRandom initialized";

    action.Run();

    LOG(INFO) << "Simulation finished.";

    // add TID tree for the generated events
    if(!cmd_noTID->isSet()) {
        LOG(INFO) << "Add TID tree to the output file";
        mc::pluto::utils::PlutoTID::AddTID(action.outfile);
    }

    return EXIT_SUCCESS;
}

void GunAction::Run() const
{
    auto nParticles = particles.size();

    // prepare particle names if verbosity level is set
    vector<string> particle_names;
    if (VLOG_IS_ON(1)) {
        std::transform(particles.begin(), particles.end(),
                       std::back_inserter(particle_names),
                       [] (const ParticleTypeDatabase::Type* p) -> std::string {
            return p->Name();
        });
    }

    VLOG(1) << "Running RandomGun";
    VLOG(1) << "number of Events:    " << nEvents;
    VLOG(1) << "number of particles: " << nParticles;
    VLOG(1) << "Particles: " << particle_names;
    VLOG(1) << "E min: " << Emin << " MeV";
    VLOG(1) << "E max: " << Emax << " MeV";
    VLOG(1) << "Theta min: "     << radian_to_degree(thetaMin)  << " degree";
    VLOG(1) << "Theta max: "     << radian_to_degree(thetaMax)  << " degree";
    if (flatTheta)
        VLOG(1) << "Generating a flat theta distribution (non-uniform x, y, and z values)";
    VLOG(1) << "opening angle: " << radian_to_degree(openAngle) << " degree";


    WrapTFileOutput file(outfile);

    TTree* tree = file.CreateInside<TTree>("data","Random Particles");

    TClonesArray* particles_array = new TClonesArray("PParticle", int(particles.size()));
    particles_array->SetOwner(kTRUE);

    tree->Branch("Particles", particles_array);

    // ALL in GeV from here!!!!

    vec3 dirFirst;
    for( unsigned evt=0; evt< nEvents; ++evt ) {

        particles_array->Clear();

        for( unsigned iParticle=0; iParticle<nParticles; ++iParticle ) {
            if (iParticle==0 || std::isnan(openAngle))
            {
                const auto randomDir = getRandomDir() * getRandomMomentum(particles.at(iParticle));
                (*particles_array)[iParticle] = new PParticle(particles.at(iParticle)->PlutoID(),randomDir);
                if (iParticle == 0)
                    dirFirst = randomDir;
            }
            else{
                const auto dirCone = getRandomDirInCone(dirFirst) * getRandomMomentum(particles.at(iParticle));
                (*particles_array)[iParticle] = new PParticle(particles.at(iParticle)->PlutoID(),dirCone);
            }
        }
        tree->Fill();
    }
}
