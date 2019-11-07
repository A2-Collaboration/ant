#include "analysis/physics/etaprime/etaprime_dalitz.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/system.h"
#include "base/ParticleType.h"

#include "analysis/plot/HistogramFactory.h"
#include "analysis/utils/ParticleTools.h"
#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "TSystem.h"
#include "TRint.h"
#include "TROOT.h"

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"


using namespace ant;
using namespace std;
using namespace RooFit;


string concat_string(const vector<string>& strings, const string& delimiter = ", ")
{
    if (strings.empty())
        return "";

    return accumulate(next(strings.begin()), strings.end(), strings.front(),
            [&delimiter] (string& concat_str, const string& str) {
                return concat_str + delimiter + str;
            });
}

string cuts_path(const vector<string>& cuts, const char* delimiter = "/")
{
    if (cuts.empty())
        return "";

    stringstream s;
    copy(cuts.begin(), cuts.end(), ostream_iterator<string>(s, delimiter));
    return s.str();
}

string get_path(const string& cut_string, const string& tree)
{
    return tree + "/" + cut_string;
}


struct q2_bin_cut_t {
    string q2_bin;
    vector<string> cuts;
    string cut_string;

    void create_cut_string() {
        cut_string = cuts_path(cuts, "/");
    }

    q2_bin_cut_t(const string& bin, const string& cuts) : q2_bin(bin), cut_string(cuts) {}

    q2_bin_cut_t(const string& bin, const vector<string>& _cuts) : q2_bin(bin), cuts(_cuts) {
        create_cut_string();
    }
};


void test_path_building()
{
    const vector<string> cuts = {
        "selection",
        "KinFitProb > 0.1",
        "nothing",
        "thight cluster size"};
    const string tree = "EtapDalitz_plot_Sig";
    cout << "Test building tree path from cuts vector:" << endl
        << get_path(concat_string(cuts, "/"), tree) << endl;
    cout << "Test cuts path using copy and stringstream:" << endl
        << cuts_path(cuts) << endl;
}

void traverseCuts(TDirectory* dir, vector<vector<string>>& cuts) {
    auto keys = dir->GetListOfKeys();
    if (!keys)
        return;

    vector<string> dirnames;
    bool h_found = false;
    TIter nextk(keys);
    TKey* key;
    TKey* nextdir = nullptr;
    while ((key = static_cast<TKey*>(nextk())))
    {
        auto classPtr = TClass::GetClass(key->GetClassName());
        if (classPtr->InheritsFrom(TDirectory::Class())) {
            const string dirname(key->GetName());
            if (dirname == "h")
                h_found = true;
            else {
                nextdir = key;
                dirnames.emplace_back(dirname);
            }
        }
    }

    if (h_found && !dirnames.empty()) {
        cuts.emplace_back(dirnames);
        if (nextdir) {
            traverseCuts(dynamic_cast<TDirectory*>(nextdir->ReadObj()), cuts);
        }
    }
}

vector<vector<string>> extractCuts(const string& prefix, const WrapTFileInput& input) {
    TDirectory* prefixDir = nullptr;
    if (!input.GetObject(prefix, prefixDir))
        throw runtime_error("Cannot find prefix dir " + prefix);
    vector<vector<string>> cuts;
    traverseCuts(prefixDir, cuts);
    return cuts;
}

void print_extracted_cuts(const string& file)
{
    WrapTFileInput input(file);
    const string prefix = "EtapDalitz_plot_Sig";
    auto cuts = extractCuts(prefix, input);
    cout << "Extracted Cuts:" << endl;
    size_t cut_level = 0;
    for (const auto& vec : cuts) {
        cout << "  cut level " << ++cut_level << endl;
        for (const auto& cut : vec)
            cout << "    " << cut << endl;
    }
}



struct TCLAPInterval : interval<int> {
    using interval::interval;
    using ValueCategory = TCLAP::ValueLike;
};

int main(int argc, char** argv) {
    SetupLogger();

    // parse command line
    TCLAP::CmdLine cmd("EtapDalitz_fit", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_debug = cmd.add<TCLAP::MultiSwitchArg>("","debug","Enable debug mode",false);

    auto cmd_ref = cmd.add<TCLAP::MultiSwitchArg>("r","reference","Run Reference Channel Analysis", false);
    auto cmd_ref_only = cmd.add<TCLAP::MultiSwitchArg>("","ref-only","Only Reference Channel Analysis", false);

    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","ROOT input file",true,"","rootfile");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");
    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup by name",true,"", &allowedsetupnames);
    auto cmd_EPTrange = cmd.add<TCLAP::ValueArg<TCLAPInterval>>("c","EPTrange","EPT channel range for reference fits, e.g. 0-40",
                                                                false,TCLAPInterval{0,40},"channels");

    cmd.parse(argc, argv);

    const bool ref = cmd_ref->isSet();
    const bool ref_only = cmd_ref_only->isSet();

    // verbosity management
    if (cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    // do some tests in the beginning to make sure all functions work as expected
    if (cmd_debug->isSet()) {
        test_path_building();
        cout << "\nCall cut extraction method\n" << endl;
        print_extracted_cuts(cmd_input->getValue());
    }

    ExpConfig::Setup::SetByName(cmd_setup->getValue());

    WrapTFileInput input(cmd_input->getValue());

    const auto taggChRange = cmd_EPTrange->getValue();
    if (!taggChRange.IsSane()) {
        LOG(ERROR) << "Provided Tagger channel range " << taggChRange << " is not sane.";
        return EXIT_FAILURE;
    }
    if (cmd_EPTrange->isSet()) {
        LOG(WARNING) << "Using non-default Tagger channel range, may not yield correct results (debugging purposes)";
    }

    // create TRint as RooFit internally creates functions/histograms,
    // prevents this stupid gStyle=0 related error, sigh...
    argc=0; // prevent TRint to parse any cmdline
    TRint app("EtapDalitz_fit",&argc,argv,nullptr,0,true);
    if (cmd_batchmode->isSet())
        gROOT->SetBatch(true);

    unique_ptr<WrapTFileOutput> masterFile;
    if (cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    if (ref || ref_only) {
        TH1* ref;
        input.GetObject("EtapDalitz_plot_Ref/KinFitProb > 0.01/PID E cut < 0.3 MeV/h/Data/etapIM_kinfitted", ref);

        // fit ARGUS model with gaus for reference
        // --- Observable ---
        RooRealVar mes("IM","IM_{#gamma#gamma} [MeV]", 840, 1020);

        // --- Parameters ---
        RooRealVar sigmean("sigmean","#eta' mass", 958., 920., 980.);
        RooRealVar sigwidth("sigwidth","#eta' width", 10., .1, 30.);

        // --- Build Gaussian PDF ---
        RooGaussian signalModel("signal","signal PDF",mes,sigmean,sigwidth);

        // --- Build Argus background PDF ---
        RooRealVar argpar("argpar","argus shape parameter",-5.,-25.,5.);
        RooArgusBG background("background","Argus PDF",mes,RooConst(1000),argpar);

        // --- Construct signal+background PDF ---
        RooRealVar nsig("nsig","#signal events",1000,0.,1e5);
        RooRealVar nbkg("nbkg","#background events",1000,0.,1e5);
        RooAddPdf model("model","g+a",RooArgList(signalModel,background),RooArgList(nsig,nbkg));

        RooDataHist h_roo_data("h_roo_data","dataset",mes,ref);

        // --- Perform extended ML fit of composite PDF to toy data ---
        model.fitTo(h_roo_data);

        // --- Plot toy data and composite PDF overlaid ---
        RooPlot* mesframe = mes.frame();
        h_roo_data.plotOn(mesframe);
        model.plotOn(mesframe);
        model.plotOn(mesframe, Components(background), LineStyle(ELineStyle::kDashed));

        mesframe->Draw();
        gPad->Modified();
        gPad->Update();
    }


    // run TRint
    if (!cmd_batchmode->isSet()) {
        if (!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {
            if (masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if (masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return EXIT_SUCCESS;
}
