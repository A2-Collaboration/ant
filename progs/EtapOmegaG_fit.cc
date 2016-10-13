#include "base/Logger.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ParticleType.h"

#include "root-addons/analysis_codes/Math.h"

#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"

#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace std;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("EtapOmegaG_fit", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    WrapTFileInput input(cmd_input->getValue());

    const vector<string> histpaths{
        "h_IM_4g_omega_PID",
        "EtapSergey/EtapSergey/Simple/Prob/IM_Omega/CBVetoSumE/gNonPi0_2/h/Data/h_IM_4g",
        "EtapOmegaG/SigPi0/DiscardedEk=0/AntiPi0FitProb<0.00001||nan/AntiEtaFitProb<0.0001||nan/TreeFitProb>0.2/gNonPi0_2/CBSumVetoE<0.2/IM_Pi0g[1]/h/Data/h_IM_4g",
    };

    TH1D* h = nullptr;

    for(auto& histpath : histpaths) {
        if(input.GetObject(histpath, h)) {
            LOG(INFO) << "Found " << histpath << " in inputfile";
            break;
        }
    }

    if(!h) {
        LOG(ERROR) << "No histogram found";
        return 1;
    }


    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    const auto npx = 500;
    const auto r_min =  870.0;
    const auto r_max = 1100.0;
    const auto etap_mass = ParticleTypeDatabase::EtaPrime.Mass();
    const auto expected_width = 5.0;


    TF1* sig = new TF1("sig", "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    // position
    sig->SetParameter(1, etap_mass);

    // width
    sig->SetParameter(2, expected_width);

    TF1* bg = new TF1("bg", "pol2", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");

    TFSum::FitRanged(h, bg, 900, 940, 980, 1100);
    bg->FixParameter(0, bg->GetParameter(0));
    bg->FixParameter(1, bg->GetParameter(1));
    bg->FixParameter(2, bg->GetParameter(2));


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
            TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

            if(masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            h->SetStats(0);
            h->Draw();
            h->Fit(sum->Function(), "REM0NB");

            sum->SyncToFcts();

            sum->Draw();

            const double peak_pos  = sig->GetParameter(1);
            const double fitted_sigma = sig->GetParameter(2);

            const auto r = interval<double>::CenterWidth(peak_pos, 6*fitted_sigma);

            const double total_area = sum->Function()->Integral(r.Start(), r.Stop());
            const double bg_area    = bg->Integral(r.Start(), r.Stop());
            const double sig_area   = total_area - bg_area;

            const double sig_to_bg = sig_area / bg_area;

            cout << "Mass offset = " << peak_pos - etap_mass << " MeV\n";
            cout << "Sig/BG      = " << sig_to_bg << "\n";
            cout << "Sig         = " << sig_area << endl;

            auto pt = new TPaveText(.1,.5,.4,.9,"NDC");
            auto add_text = [pt] (const string& txt) {
                pt->AddText(txt.c_str());
            };
            add_text("N="+to_string(sig_area));
            add_text("#sigma="+to_string(fitted_sigma));
            add_text("#mu="+to_string(peak_pos));
            pt->Draw();

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return 0;
}
