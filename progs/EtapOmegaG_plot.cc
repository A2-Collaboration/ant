#include "base/Logger.h"

#include "TRint.h"
#include "analysis/plot/root_draw.h"
#include "analysis/physics/etaprime/etaprime_omega_gamma.h"

#include "base/CmdLine.h"
#include "base/interval.h"

#include "base/printable.h"
#include "base/iterators.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"

#include "TH1D.h"

#include <string>
#include <list>

using namespace ant;
using namespace std;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_file = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");;
    cmd.parse(argc, argv);

    WrapTFileInput input(cmd_file->getValue());

    TTree* t;
    if(!input.GetObject("EtapOmegaG/treeCommon",t))
        return 1;

    analysis::physics::EtapOmegaG::TreeCommon treeCommon;
    treeCommon.LinkBranches(t);

    if(!input.GetObject("EtapOmegaG/treeSig",t))
        return 1;
    analysis::physics::EtapOmegaG::Sig_t::Tree_t treeSig;
    treeSig.LinkBranches(t);

    for(int i=0;i<1000;i++) {
        treeSig.Tree->GetEntry(i);
        treeCommon.Tree->GetEntry(i);

        cout << treeCommon.CBSumE << endl;
        cout << treeSig.ggg() << endl;
    }

    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("plot",&fake_argc,fake_argv,nullptr,0,true);

    app.Run();

    return 0;
}