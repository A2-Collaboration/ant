#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"

#include "tree/TAntHeader.h"

#include "PParticle.h"
#include "PReaction.h"

#include "TSystem.h"
#include "TRint.h"

#include <list>

#include "base/vec/LorentzVec.h"
#include "base/vec/vec3.h"
#include "base/ParticleType.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "base/std_ext/math.h"

using namespace ant;

using namespace std;

volatile static bool interrupt = false;

struct ThetaPhi_t {
    double theta;
    double phi;
    ThetaPhi_t(const double t=0.0, const double p=0.0): theta(t), phi(p) {}
};

ThetaPhi_t getRandomThetaPhi() {
    double rnd[2];
    gRandom->RndmArray(2,rnd);
    const double theta = TMath::ACos(rnd[0]);
    const double phi   = 2.0 * TMath::Pi() * rnd[1];
    return {theta,phi};
}

std::pair<LorentzVec,LorentzVec> decayIsotropicallyCMS (const double m) {

    const auto tp = getRandomThetaPhi();
    const auto   p = vec3::RThetaPhi(m/2.0, tp.theta, tp.phi);
    return {{p, m/2},{-p,m/2}};
}



vec3 getRandomDir()
{
    vec3 dir;
    gRandom->Sphere(dir.x, dir.y, dir.z,1.0);
    return dir;
}


double getRandomMomentum(const double mass, const interval<double> Erange)
{
    const auto m = mass;
    const auto E = gRandom->Uniform(Erange.Start(), Erange.Stop()) + m;
    const auto p = sqrt(E*E - m*m);
    return p;
}

LorentzVec rndm(const double mass, const interval<double> Erange)
{
    const auto m = mass;
    const auto E = gRandom->Uniform(Erange.Start(), Erange.Stop()) + m;
    const auto p = sqrt(E*E - m*m);
    return {{p*getRandomDir()}, E};
}

// energies in GeV
bool similar(const double a, const double b) {
    return int(a*1000.0) / 25 == int(b*1000.0) / 25;
}

int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true;} );

    constexpr auto GeV = 1000.0;

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_output   = cmd.add<TCLAP::ValueArg<string>>("o","output",   "Output file",      true,"","filename");
//    auto cmd_particle = cmd.add<TCLAP::ValueArg<string>>("p","particle", "Particle to decay",false,"pi0","particle");
//    auto cmd_reaction = cmd.add<TCLAP::ValueArg<string>>("r","reaction", "Reaction string",  false,"g g","reaction");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_Emin      = cmd.add<TCLAP::ValueArg<double>>      ("",  "Emin",         "Minimal incident energy [MeV]", false, 0.0, "double [MeV]");
    auto cmd_Emax      = cmd.add<TCLAP::ValueArg<double>>      ("",  "Emax",         "Maximal incident energy [MeV]", false, 1.6*GeV, "double [MeV]");
    auto cmd_events    = cmd.add<TCLAP::ValueArg<int>>         ("n",  "",            "number of events", false, 10000, "n");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    const auto Erange = interval<double>(cmd_Emin->getValue(), cmd_Emax->getValue()) / 1000.0;

    WrapTFileOutput outfile(cmd_output->getValue(), true);

    TTree* tree = new TTree("Particles","");

    const auto mass = ParticleTypeDatabase::Pi0.Mass() / 1000.0; // GeV

    TClonesArray* storage = new TClonesArray("PParticle", 3);

    tree->Branch("data", storage);

    PParticle* pi0 = new PParticle("pi0");
    PParticle* g1  = new PParticle("g");
    PParticle* g2  = new PParticle("g");
    (*storage)[0] = pi0;
    (*storage)[1] = g1;
    (*storage)[2] = g2;

    g1->SetParentId(pi0->ID());
    g2->SetParentId(pi0->ID());

    g1->SetParentIndex(0);
    g2->SetParentIndex(0);

    const int nevents = cmd_events->getValue();
    while(tree->GetEntries() < nevents) {
        auto photons = decayIsotropicallyCMS(mass);

        const auto p = getRandomMomentum(mass, Erange) * getRandomDir();
        const auto pi0lv = rndm(mass,Erange);
        {
            const auto boost = pi0lv.BoostVector();
            photons.first.Boost(boost);
            photons.second.Boost(boost);
        }


        if(similar(photons.first.E, photons.second.E)) {

            pi0->SetVect4(pi0lv);

            g1->SetVect4(photons.first);
            g2->SetVect4(photons.second);

            tree->Fill();

        }


    }



    return EXIT_SUCCESS;
}
