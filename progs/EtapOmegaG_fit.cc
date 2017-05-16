#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/ParticleType.h"

#include "analysis/plot/RootDraw.h"
#include "root-addons/analysis_codes/Math.h"

#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"

#include "TSystem.h"
#include "TRint.h"
#include "Math/Interpolator.h"

using namespace ant;
using namespace std;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("EtapOmegaG_fit", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_data = cmd.add<TCLAP::ValueArg<string>>("","data","Data input",true,"","rootfile");
    auto cmd_mc = cmd.add<TCLAP::ValueArg<string>>("","mc","MC signal/reference input",true,"","rootfile");
    auto cmd_histpath = cmd.add<TCLAP::ValueArg<string>>("","histpath","Path for hists (determines cutstr)",false,"EtapOmegaG_plot_Ref/DiscardedEk=0/KinFitProb>0.02","path");
    auto cmd_histname = cmd.add<TCLAP::ValueArg<string>>("","histname","Name of hist",false,"h_IM_2g","name");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    TH1D* h_data = nullptr;
    WrapTFileInput input_data(cmd_data->getValue()); // keep it open
    {
        const string histpath = cmd_histpath->getValue()+"/h/Data/"+cmd_histname->getValue();
        if(!input_data.GetObject(histpath, h_data)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    TH1D* h_mc = nullptr;
    WrapTFileInput input_mc(cmd_mc->getValue()); // keep it open
    {
        const string histpath = cmd_histpath->getValue()+"/h/Sum_MC/"+cmd_histname->getValue();
        if(!input_mc.GetObject(histpath, h_mc)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    const auto npx = 500;
    const unsigned polorder = 6;
    const double xmin = h_data->GetXaxis()->GetXmin();
    const double xmax = h_data->GetXaxis()->GetXmax();


    struct sig_func_t {

        // interpolator can't be copied, so create it on heap
        std::shared_ptr<ROOT::Math::Interpolator> Interpolator;
        const double xmin;
        const double xmax;
//        double x_maxpos = 0.0;

        explicit sig_func_t(const TH1D& h) :
            Interpolator(make_shared<ROOT::Math::Interpolator>()),
            xmin(h.GetXaxis()->GetBinCenter(1)),
            xmax(h.GetXaxis()->GetBinCenter(h.GetNbinsX()))
        {
            vector<double> x;
            vector<double> y;
            const double integral = h.Integral(1, h.GetNbinsX());
            for(int bin=1;bin<=h.GetNbinsX();bin++) {
                x.push_back(h.GetXaxis()->GetBinCenter(bin));
                y.push_back(h.GetBinContent(bin)/integral);
            }
            Interpolator->SetData(x, y);

//            TF1 f_interp("temp",[this] (double* x, double*) {
//                return Interpolator->Eval(x[0]);
//            },xmin,xmax,0);
//            f_interp.SetNpx(1000);
//            x_maxpos = f_interp.GetMaximumX(xmin, xmax);
        }

        double operator()(double* x_, double* p) const {
            const double x  = x_[0];
            const double N  = p[0];
            const double x0 = p[1];
            const double s  = p[2];
            const double x_shiftscale = s*(x-x0)+x0;
            if(x_shiftscale<xmin)
                return 0;
            if(x_shiftscale>xmax)
                return 0;
            return N*Interpolator->Eval(x_shiftscale);
        }
    };

    sig_func_t sig_func(*h_mc);

    auto sig = new TF1("sig", sig_func, xmin, xmax, 3);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // normalization
    sig->SetParameter(0, h_mc->Integral(xmin, xmax));
    sig->SetParName(0, "N");

    // position shift
    sig->SetParameter(1, 958);
    sig->SetParName(1, "x_0");

    // width scale
    sig->SetParameter(2, 1.0);
    sig->SetParName(2, "s");

//    sig->FixParameter(2, sig->GetParameter(2));

    auto bg = new TF1("bg", ("pol"+to_string(polorder)).c_str(), xmin, xmax);
    bg->SetLineColor(kBlue);
    for(unsigned i=0;i<=polorder;i++) {
        bg->SetParameter(i,0.0);
        bg->SetParName(i, ("BG p"+to_string(i)).c_str());
    }

    TFSum::FitRanged(h_data, bg, xmin, 930, 980, xmax);

    TFSum sum("sum", sig, bg, xmin, xmax);
    sum.SetNpx(npx);

    h_data->Fit(sum.Function(), "REM0NB");
//    h_mc->Fit(sig, "REM0NB");
    sum.SyncToFcts();


    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
            TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

            if(masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            canvas("EtapOmegaG_fit")
                    << h_data << samepad << sum
//                    << h_mc << samepad << sig
                    << endc;

//            new TCanvas();
//            sig->Draw();


//            h->SetStats(0);
//            h->Draw();
//            h->Fit(sum->Function(), "REM0NB");

//            sum->SyncToFcts();

//            sum->Draw();

//            const double peak_pos  = sig->GetParameter(1);
//            const double fitted_sigma = sig->GetParameter(2);

//            const auto r = interval<double>::CenterWidth(peak_pos, 6*fitted_sigma);

//            const double total_area = sum->Function()->Integral(r.Start(), r.Stop());
//            const double bg_area    = bg->Integral(r.Start(), r.Stop());
//            const double sig_area   = total_area - bg_area;

//            const double total_area_r = sum->Function()->Integral(peak_pos, r.Stop());
//            const double bg_area_r    = bg->Integral(peak_pos, r.Stop());
//            const double sig_area_r   = total_area_r - bg_area_r;

//            const double sig_to_bg = sig_area / bg_area;

//            cout << "Mass offset = " << peak_pos - etap_mass << " MeV\n";
//            cout << "Sig/BG      = " << sig_to_bg << "\n";
//            cout << "Sig         = " << sig_area << endl;
//            cout << "Sig_r       = " << sig_area_r << endl;


//            auto pt = new TPaveText(.1,.5,.4,.9,"NDC");
//            auto add_text = [pt] (const string& txt) {
//                pt->AddText(txt.c_str());
//            };
//            add_text("N="+to_string(sig_area));
//            add_text("N_r="+to_string(sig_area_r));
//            add_text("#sigma="+to_string(fitted_sigma));
//            add_text("#mu="+to_string(peak_pos));
//            pt->Draw();

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return 0;
}
