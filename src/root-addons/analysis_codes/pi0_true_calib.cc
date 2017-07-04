#include "pi0_true_calib.h"
#include "base/Logger.h"
#include "analysis/plot/HistogramFactory.h"

#include "TH1D.h"
#include "analysis/plot/RootDraw.h"
#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"
#include "TCanvas.h"
#include "tree/TAntHeader.h"
#include <TMultiGraph.h>
#include "PParticle.h"
#include "PReaction.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TRint.h"
#include "TRandom3.h"
#include <list>
#include "base/vec/LorentzVec.h"
#include "base/vec/vec3.h"
#include "base/ParticleType.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "base/std_ext/math.h"
#include <math.h>

TRandom3 Random;

using namespace ant;
using namespace std;
using namespace ant::analysis;

volatile static bool interrupt = false;

static interval<double> Erange = {0,0};

struct ThetaPhi_t {
    double theta;
    double phi;
    ThetaPhi_t(const double t=0.0, const double p=0.0): theta(t), phi(p) {}
};

ThetaPhi_t getRandomThetaPhi() {
    double rnd[2];
    Random.RndmArray(2,rnd);
    const double theta = TMath::ACos(rnd[0]);
    const double phi   = 2.0 * TMath::Pi() * rnd[1];
    return {theta,phi};
}

ThetaPhi_t getzBoostThetaPhi()
{
    auto tp = getRandomThetaPhi();
    while(tp.theta > 0.34) //20 Degree
    {
        tp = getRandomThetaPhi();
    }
    return tp;
}

std::pair<LorentzVec,LorentzVec> decayIsotropicallyCMS (const double m) {

    const auto tp = getRandomThetaPhi();
    const auto   p = vec3::RThetaPhi(m/2.0, tp.theta, tp.phi);
    return {{p, m/2},{-p,m/2}};
}


vec3 getRandomDir()
{

    vec3 dir;
    Random.Sphere(dir.x, dir.y, dir.z,1.0);
    return dir;
}

vec3 getzboost()
{
    const auto tp = getzBoostThetaPhi();

    vec3 dir = vec3::RThetaPhi(1.0, tp.theta, tp.phi);

    return dir;
}


LorentzVec rndPhotonbeamenergy()
{
    LorentzVec beam;// = {{x,y,z},E};
    const auto E=Random.Uniform(Erange.Start(), Erange.Stop());
    beam.E = E;
    beam.p = vec3(0,0,E);
    return beam;
}


std::pair<LorentzVec,LorentzVec> Pi0Boost()
{

    auto m_pi = ParticleTypeDatabase::Pi0.Mass()/1000.0;
    auto m_p  = ParticleTypeDatabase::Proton.Mass()/1000.0;

    const LorentzVec Target = {{0,0,0}, m_p};
    const LorentzVec beam = rndPhotonbeamenergy();

    const LorentzVec BT = Target + beam;
    const auto invmass = BT.M();

    LorentzVec Pi0;
    LorentzVec Proton;

    auto p    = sqrt(((pow(m_pi,4) - 2 * pow(m_pi,2) * pow(m_p,2) - 2 * pow(m_pi,2) * pow((invmass),2) + pow(m_p,4) - 2 * pow(m_p , 2) * pow((invmass),2) + pow((invmass),4))/(4 *pow((invmass),2))));
    Pi0.E     = sqrt(pow(m_pi,2)+ pow(p,2));
    Proton.E  = sqrt(pow(m_p,2)+ pow(p,2));


    vec3 tp = getRandomDir();
    Pi0.p = vec3::RThetaPhi(p,tp.Theta(),tp.Phi());
    Proton.p = - Pi0.p;

    Pi0.Boost(BT.BoostVector());

    return {Pi0,Proton};
}




double getRandomMomentum(const double mass, const interval<double> Erange)
{
    const auto m = mass;
    const auto E = Random.Uniform(Erange.Start(), Erange.Stop()) + m;
    const auto p = sqrt(E*E - m*m);
    return p;
}

LorentzVec rndm(const double mass, bool opt)
{

    const auto E = Random.Uniform(Erange.Start(), Erange.Stop()) + mass;
    const auto p = sqrt(E * E - mass*mass);
    if(opt==false){
    return {{p*getRandomDir()}, E};

    }
    else {
        return {{p*getzboost()}, E};
    }
}

// energies in GeV
bool similar(const double a, const double b) {
    double ediff = a - b;
    if(ediff < 0)
    {
        ediff *= -1;
    }

    return ediff<=0.025;
}

void pi0_true_calib::Do()
{


    HistogramFactory HistFac("h");

    const BinSettings bins_IM   (400, 0, 1100); // MeV


    auto h_IM_CB  = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800));
    auto h_IM_Shifted_Target  = HistFac.makeTH2D("Decay not in the origin",   "Target [cm]","E_{#gamma} [MeV]",BinSettings(100,-5,5),bins_IM);






//    LOG(INFO) << "Hej!";

    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true;} );

//    constexpr auto GeV = 1000.0;

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
//    auto cmd_output   = cmd.add<TCLAP::ValueArg<string>>("o","output",   "Output file",      true,"","filename");
//    auto cmd_particle = cmd.add<TCLAP::ValueArg<string>>("p","particle", "Particle to decay",false,"pi0","particle");
//    auto cmd_reaction = cmd.add<TCLAP::ValueArg<string>>("r","reaction", "Reaction string",  false,"g g","reaction");
//    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
//    auto cmd_Emin      = cmd.add<TCLAP::ValueArg<double>>      ("",  "Emin",         "Minimal incident energy [MeV]", false, 0.0, "double [MeV]");
//    auto cmd_Emax      = cmd.add<TCLAP::ValueArg<double>>      ("",  "Emax",         "Maximal incident energy [MeV]", false, 1.6*GeV, "double [MeV]");
//    auto cmd_events    = cmd.add<TCLAP::ValueArg<int>>         ("n",  "",            "number of events", false, 10000, "n");
//    auto cmd_reqsym    = cmd.add<TCLAP::SwitchArg>             ("",   "sym",         "Require symmetric photon energies");
//    auto cmd_zboost    = cmd.add<TCLAP::SwitchArg>             ("",   "zboost",      "Boost the Pions in z-Direction; True or False");
//    auto cmd_Prod      = cmd.add<TCLAP::SwitchArg>             ("",   "Prod",        "Get the Product of the Pion; Change Beam Energy with E_min and E_max"  );


    Random.SetSeed();

//    cmd.parse(argc, argv);
//    if(cmd_verbose->isSet()) {
//        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
//    }

    const bool sym = true;
    bool zboost = true;
    bool Prod =true;

    Erange = interval<double>(1480, 1604) / 1000.0;

//    WrapTFileOutput outfile(cmd_output->getValue(), true);

    TTree* tree = new TTree("data","");

    const auto mass = ParticleTypeDatabase::Pi0.Mass() / 1000.0; // GeV

/*    const auto nParticles = Prod ? 4 : 3*/;

    TClonesArray* storage = new TClonesArray("PParticle", 3);

    tree->Branch("Particles", storage);

    PParticle* pi0 = new PParticle("pi0");
    PParticle* g1  = new PParticle("g");
    PParticle* g2  = new PParticle("g");
    (*storage)[0] = pi0;
    (*storage)[1] = g1;
    (*storage)[2] = g2;

    g1->SetParentId(pi0->ID());
    g2->SetParentId(pi0->ID());
    pi0->SetDaughterIndex(1);

    g1->SetParentIndex(0);
    g2->SetParentIndex(0);

//    PParticle* proton = nullptr;
//    if(Prod) {
//        proton = new PParticle("p");
//        (*storage)[3] = proton;
//        proton->
//    }


    const int nevents = 100000;
    while(tree->GetEntries() < nevents) {

        LorentzVec pi0lv;

        if(Prod==false){

            pi0lv = rndm(mass,zboost);
        }
        else
        {
            const auto pip = Pi0Boost();
            pi0lv = pip.first;
        }

        auto photons = decayIsotropicallyCMS(mass);
        {
            const auto boost = pi0lv.BoostVector();
            photons.first.Boost(boost);
            photons.second.Boost(boost);
        }


        if(!sym || similar(photons.first.E, photons.second.E)) {

            pi0->SetVect4(pi0lv);

            g1->SetVect4(photons.first);
            g2->SetVect4(photons.second);

            auto CB_shower_radius = 40;
            auto targetshift = Random.Uniform(-5,5);



            auto theta_true = g1->Vect().Angle(g2->Vect());
            auto m_pi0_true = sqrt(2 * g1->E()*g2->E() * (1-cos(theta_true))) * 1000;


            auto theta_shifted = atan(sin(theta_true)*CB_shower_radius / (cos(theta_true) * (CB_shower_radius + targetshift)));
            auto m_pi0_shift = sqrt(2 * g1->E()*g2->E() * (1-cos(theta_shifted))) * 1000;


            tree->Fill();
            h_IM_CB->Fill(m_pi0_true,g1->E() * 1000);
            h_IM_CB->Fill(m_pi0_true,g2->E() * 1000);

            h_IM_Shifted_Target->Fill(targetshift,m_pi0_shift);
            h_IM_Shifted_Target->Fill(targetshift,m_pi0_shift);

        }


    }
    canvas()<<h_IM_CB<<h_IM_Shifted_Target<<endc ;


}
