#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "TROOT.h"
#include "TRint.h"
#include <iostream>
#include <cstring>
#include "TTree.h"
#include "TClonesArray.h"
#include "PParticle.h"
#include "analysis/utils/MCWeighting.h"
#include "analysis/plot/HistogramFactory.h"
using namespace std;
using namespace ant;
using namespace ant::calibration::gui;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-mcsmearing", ' ', "0.1");
    auto cmd_verbose  = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_input    = cmd.add<TCLAP::ValueArg<string>>("i","input",  "Input pluto file", true, "", "pluto.root");
    auto cmd_output   = cmd.add<TCLAP::ValueArg<string>>("o","output",  "Output root file",true, "","");
    auto cmd_meson    = cmd.add<TCLAP::ValueArg<string>>("m","meson",  "Meson ", true, "", "omega");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    WrapTFileInput infile(cmd_input->getValue());

    TTree* intree = nullptr;
    infile.GetObject("data",intree);
    TClonesArray* buffer = nullptr;

    intree->SetBranchAddress("Particles", &buffer);

    WrapTFileOutput outfile(cmd_output->getValue());
    TTree* outtree = outfile.CreateInside<TTree>("data", "");
    outtree->Branch("Particles", buffer);

    const analysis::utils::MCWeighting::item_t* mcitem = nullptr;
    double norm = 1.0;
    int pluto_pid=0;
    if(cmd_meson->getValue() == "omega") {
        mcitem =  &analysis::utils::MCWeighting::Omega;
        norm = 6.91090694237991715e-02;
        pluto_pid = 52;
    } else
    if(cmd_meson->getValue() == "pi0") {
        mcitem =  &analysis::utils::MCWeighting::Pi0;
        norm =4.90141090503409238e+00;
        pluto_pid = 7;
    } else {
        LOG(FATAL) << "Wromg Meson";
    }

    analysis::HistogramFactory hf("mcw");
    analysis::utils::MCWeighting mcw(hf, *mcitem);

    const auto findOneParticle = [] (const TClonesArray& c, const int pid) {
        for(int i=0;i<c.GetEntries(); ++i) {
            PParticle* p = static_cast<PParticle*>(c.At(i));
            if(p->ID() == pid) {
                return p;
            }
        }
        return (PParticle*)nullptr;
    };

    TRandom3 rng(0);
    Long64_t accepted = 0;

    for(Long64_t i=0; i<intree->GetEntries(); ++i) {
        intree->GetEntry(i);

        const auto p = findOneParticle(*buffer, pluto_pid);
        const auto bp = findOneParticle(*buffer, 14001);
        if(p && bp) {
            const auto E = bp->E()*1000.0 - ParticleTypeDatabase::Proton.Mass(); //MeV
            LorentzVec boosted = *p;
            boosted.Boost(-bp->BoostVector());
            const auto cost = cos(boosted.Theta());
            const auto rw = mcw.GetN(E, cost) / norm;
            if(rw > 1.0 || rng.Uniform(0.0,1.0) <= rw) {
                outtree->Fill();
                ++accepted;
            }
        }

    }

    LOG(INFO) << "In:" << intree->GetEntries() << ", out: " << accepted << " = " << accepted / double(intree->GetEntries()) * 100.0 << " %";

    return EXIT_SUCCESS;
}
