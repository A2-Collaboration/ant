/**
  * @file PCA.cc
  * @brief Tool to perform a Principal Component Analysis, based on ROOTs TPrincipal class, to create principal components of a given data set
  */


#include "analysis/physics/etaprime/etaprime_dalitz.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"
#include "base/WrapTTree.h"

#include <vector>
#include <algorithm>
#include <assert.h>
#include <typeinfo>
#include <functional>

#include "TPrincipal.h"
#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace std;

volatile sig_atomic_t interrupt = false;


template <typename Tree_t>
struct Data_t {

    struct GetVal_t {
        const Tree_t& Tree;

        GetVal_t(const Tree_t& t) : Tree(t) {}

        double TaggW() const {
            return Tree.TaggW;
        }
    };

    // function to store lambdas which return the desired variable value
    template <typename Variable_t>
    using get_variable_t = std::function<const Variable_t(const GetVal_t&)>;

    template <typename Variable_t>
    struct VariableGetter_t {
        get_variable_t<Variable_t> getval;
        Variable_t val;
        const string name;

        VariableGetter_t(get_variable_t<Variable_t> v,
                         const string& name,
                         Variable_t value = 0):
            getval(v),
            val(value),
            name(name)
        {}

        void GetVal(const GetVal_t& data) {
            VLOG(7) << "data type: " << abi::__cxa_demangle(typeid(getval(data)).name(), nullptr, nullptr, nullptr);
            val = getval(data);
        }
    };

    template <typename Variable_t>
    struct VariableManager_t : std::vector<VariableGetter_t<Variable_t>> {
        using vector<VariableGetter_t<Variable_t>>::vector;

        std::vector<Variable_t> values;

        void GetVal(const GetVal_t& data) {
            for (auto& v : *this)
                v.GetVal(data);
        }

        const std::vector<std::string> GetNames() {
            vector<string> names;
            for (auto& v : *this)
                names.emplace_back(v.name);
            return move(names);
        }

        const std::vector<Variable_t>& data() {
            values.clear();
            for (auto& v : *this)
                values.emplace_back(v.val);
            return cref(values);
        }
    };

    VariableManager_t<double> v;

    Data_t() {}

    void AddVariable(const string name, get_variable_t<double> g) {
        v.emplace_back(VariableGetter_t<double>(g, name));
    }

    void GetVal(const GetVal_t& data) {
        v.GetVal(data);
    }

    size_t get_number_variables() const {
        return v.size();
    }

    const vector<string> get_names() {
        return v.GetNames();
    }

    void print() {
        auto names = v.GetNames();
        auto vals = v.data();
        VLOG(5) << "current variable values:";
        for (size_t i = 0; i < v.size(); i++)
            VLOG(5) << " " << names.at(i) << ": " << vals.at(i);
    }

    const double* data() {
        return &(v.data())[0];
    }
};

struct SigData_t : Data_t<physics::EtapDalitz::SigTree_t> {
    using Tree_t = physics::EtapDalitz::SigTree_t;
    using GetVal_t = Data_t<Tree_t>::GetVal_t;


    SigData_t() : Data_t() {
        AddVariable("kinfitted eta' mass", [] (const GetVal_t& g) {
            return g.Tree.etap_kinfit().M();
        });
    }

};


int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Principal Component Analysis", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i", "input", "Input file", true, "", "input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b", "batch", "Run in batch mode (no ROOT shell afterwards)", false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m", "maxevents", "Process only max events", false, "maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o", "output", "Output file", false, "", "filename");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    cmd.parse(argc, argv);

    if (cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());


    // open the TRint app as early as possible to prevent ROOT to create a new one automatically
    // which will cause problems because of bad ROOT internal data / pointer handling, might cause segfaults
    argc = 0;  // prevent TRint to parse any cmdline
    TRint app("PCA", &argc, argv, nullptr, 0, true);


    // set signal handler after starting TRint, otherwise it will be overwritten with ROOT handlers
    signal(SIGINT, [] (int) {
        LOG(WARNING) << "Processing interrupted";
        interrupt = true;
    });


    WrapTFileInput input;
    //string setup_name;
    try {

        input.OpenFile(cmd_input->getValue());

//        auto header = input.GetSharedClone<TAntHeader>("AntHeader");

//        if (!header) {
//            LOG(WARNING) << "No TAntHeader found in " << cmd_input->getValue();
//            return 1;
//        }

//        setup_name = header->SetupName;

    } catch (const std::runtime_error& e) {
        LOG(ERROR) << "Can't open " << cmd_input->getValue() << " " << e.what();
    }


    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if (!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree " + treename + " in input file");
        if (expected_entries >= 0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree " + treename + " does not have entries == " + to_string(expected_entries));
        if (wraptree->Matches(t, false)) {
            wraptree->LinkBranches(t);
            return true;
        }
        return false;
    };


    SigData_t::Tree_t sigTree;

    if (!link_branches("EtapDalitz/signal", addressof(sigTree), -1)) {
        LOG(ERROR) << "Cannot link branches of signal tree";
        return 1;
    }

    unique_ptr<WrapTFileOutput> masterFile;
    if (cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                           true,  // cd into masterFile upon creation
                                                           WrapTFileOutput::mode_t::recreate);
    }

    HistogramFactory HistFac("PCA");


    double secs_used_signal = 0.;

    {  // main part
        const auto sigEntries = sigTree.Tree->GetEntries();

        auto signal_variables = SigData_t();

        const int n = signal_variables.get_number_variables();
        VLOG(1) << n << " variables will be used";
        auto names = signal_variables.get_names();
        assert(static_cast<size_t>(n) == names.size());
        LOG(INFO) << "The following variables have been defined:";
        for (const auto& name : names)
            LOG(INFO) << " * " << name;

        TPrincipal* principal = new TPrincipal(n, "ND");

        LOG(INFO) << "Signal Tree entries = " << sigEntries;

        auto max_entries = sigEntries;
        if (cmd_maxevents->isSet() && cmd_maxevents->getValue().back() < sigEntries) {
            max_entries = cmd_maxevents->getValue().back();
            LOG(INFO) << "Running until " << max_entries;
        }

        long long entry = 0;
        double last_percent = 0;
        ProgressCounter::Interval = 3;
        ProgressCounter progress(
                    [&entry, &sigEntries, &last_percent] (std::chrono::duration<double> elapsed) {
            const double percent = 100.*entry/sigEntries;
            const double speed = (percent - last_percent)/elapsed.count();
            LOG(INFO) << setw(2) << setprecision(4) << "Processed " << percent << " %, ETA: " << ProgressCounter::TimeToStr((100-percent)/speed);
            last_percent = percent;
        });

        while (entry < max_entries) {
            if (interrupt)
                break;

            sigTree.Tree->GetEntry(entry++);
            signal_variables.GetVal({sigTree});
            auto data = signal_variables.data();
            if (el::Loggers::verboseLevel() >= 5)
                signal_variables.print();
            principal->AddRow(data);

            ProgressCounter::Tick();
        }

        LOG(INFO) << "Finished processing the tree";
        secs_used_signal = progress.GetTotalSecs();
        LOG(INFO) << "Analysed " << entry << " events, speed "
                  << entry/secs_used_signal << " event/s";

        LOG(INFO) << "Start performing the PCA";

        // Do the actual analysis
        principal->MakePrincipals();
        // Print out the result on
        principal->Print();
        // Test the PCA
        principal->Test();
        // Make some histograms of the orginal, principal, residue, etc data
        principal->MakeHistograms();
        // Make two functions to map between feature and pattern space
        principal->MakeCode();

        secs_used_signal = progress.GetTotalSecs();
        LOG(INFO) << "Finished PCA, principal components created";
        LOG(INFO) << "Total time used: " << ProgressCounter::TimeToStr(secs_used_signal);
    }

    if (!cmd_batchmode->isSet()) {
        if (!std_ext::system::isInteractive())
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        else {
            if (masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if (masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return EXIT_SUCCESS;
}
