#include <iostream>
#include <string>
#include <algorithm>

#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "base/WrapTFile.h"
#include "TROOT.h"
#include "TRint.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "PParticle.h"

#include "detail/Corrections.h"

using namespace std;
using namespace ant;
using namespace ant::progs::corrections;


int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("Ant tool to apply radiative corrections to a Pluto file via rejection", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose", "Verbosity level (0..9)", false, 0, "level");
    auto cmd_input   = cmd.add<TCLAP::ValueArg<string>>("i","input", "Input Pluto file", true, "", "pluto.root");
    auto cmd_output  = cmd.add<TCLAP::ValueArg<string>>("o","output", "Output root file", true, "", "output.root");
    auto cmd_nEvents = cmd.add<TCLAP::ValueArg<int>>("m","maxEvents", "Max. number of events to be read", false, 0, "#events");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    const int maxEvents = cmd_nEvents->getValue();
    if (maxEvents)
        LOG(INFO) << "Max. number of events set to " << maxEvents;

    WrapTFileInput input(cmd_input->getValue());

    TTree* input_tree = nullptr;
    if (!input.GetObject("data", input_tree))
        LOG(FATAL) << "\"data\" not found, make sure the provided file is a Pluto file";

    if (!input_tree->GetEntries())
        LOG(FATAL) << "The provided tree seems to be empty";

    TClonesArray* buffer = nullptr;

    int meson_id = -1, lm_id = -1, lp_id = -1;
    Channel channel = Channel::unknown;
    {
        bool dilepton = false, dimuon = false, pi0 = false, eta = false, etap = false;
        PStaticData* pluto_db = makeStaticData();
        input_tree->SetBranchAddress("Particles", &buffer);
        input_tree->GetEntry(0);

        for (int i=0; i < buffer->GetEntries(); ++i) {
            PParticle* p = static_cast<PParticle*>(buffer->At(i));
            if (p->ID() == pluto_db->GetParticleID("pi0")) {
                pi0 = true;
                meson_id = pluto_db->GetParticleID("pi0");
                channel = Channel::pi0_eeg;
            } else if (p->ID() == pluto_db->GetParticleID("eta")) {
                eta = true;
                meson_id = pluto_db->GetParticleID("eta");
            } else if (p->ID() == pluto_db->GetParticleID("eta'")) {
                etap = true;
                meson_id = pluto_db->GetParticleID("eta'");
            } else if (p->ID() == pluto_db->GetParticleID("dilepton")) {
                dilepton = true;
                lm_id = pluto_db->GetParticleID("e-");
                lp_id = pluto_db->GetParticleID("e+");
            } else if (p->ID() == pluto_db->GetParticleID("dimuon")) {
                dimuon = true;
                lm_id = pluto_db->GetParticleID("mu-");
                lp_id = pluto_db->GetParticleID("mu+");
            }
        }

        if (!pi0 && !eta && !etap)
            LOG(FATAL) << "No pseudo-scalar meson found";

        if (pi0 == eta ? pi0 : etap)
            LOG(FATAL) << "More than one pseudo-scalar meson found";

        if (!dilepton && !dimuon)
            LOG(FATAL) << "No dilepton pair found";

        if (pi0 && dimuon)
            LOG(FATAL) << "Found pi0 and mu+ mu- pair in the sample";

        if (eta)
            channel = dilepton ? Channel::eta_eeg : Channel::eta_mumug;
        else if (etap)
            channel = dilepton ? Channel::etap_eeg : Channel::etap_mumug;
    }
    Corrections corr(channel);
    const double weight_max = 1. + corr.get_max_correction()/100.;
    VLOG(1) << "max correction = " << corr.get_max_correction() << "; max weight = " << weight_max;

    input_tree->SetBranchAddress("Particles", &buffer);

    WrapTFileOutput output(cmd_output->getValue());
    TTree* output_tree = output.CreateInside<TTree>("data", "");
    output_tree->Branch("Particles", buffer);

    const auto find_particle = [] (const TClonesArray& c, const int pid) {
        for (int i=0; i < c.GetEntries(); ++i) {
            PParticle* p = static_cast<PParticle*>(c.At(i));
            if (p->ID() == pid)
                return p;
        }
        return static_cast<PParticle*>(nullptr);
    };

    TRandom3 rng(0);
    long long accepted = 0;
    long long read = 0;

    for (long long i = 0; i < input_tree->GetEntries(); ++i) {
        if (maxEvents && i >= maxEvents)
            break;

        input_tree->GetEntry(i);
        ++read;

        const auto meson = find_particle(*buffer, meson_id);
        const auto lm = find_particle(*buffer, lm_id);
        const auto lp = find_particle(*buffer, lp_id);

        if (meson && lm && lp) {
            double weight = 1. + corr.interpolate(lm, lp, meson)/100.;
            // weight normalized to weight_max, hence it's always between 0 and 1
            double w = weight/weight_max;

            // rejection
            if (rng.Uniform(0.,1.) <= w) {
                output_tree->Fill();
                ++accepted;
            }
        } else
            LOG(ERROR) << "Not all needed particles found in event " << i;
    }

    double percent_accepted = double(accepted) / read * 100.;
    LOG(INFO) << "in: " << read << ", out: " << accepted << "; "
              << percent_accepted << "% accepted, " << 100.-percent_accepted << "% rejected";

    return EXIT_SUCCESS;
}
