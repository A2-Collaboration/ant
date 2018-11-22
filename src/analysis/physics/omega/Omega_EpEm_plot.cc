#include "Omega_EpEm.h"
#include "physics/Plotter.h"

#include "analysis/plot/CutTree.h"

#include "base/Logger.h"
#include "base/interval.h"
#include "base/std_ext/vector.h"

#include "expconfig/ExpConfig.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

struct Omega_EpEm_plot : Plotter {



    Omega_EpEm_plot(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        Plotter(name, input, opts)
    {


    }

//    static void init_tree(const WrapTFileInput& input, WrapTTree& tree, const string& name)
//    {

//    }

//    virtual long long GetNumEntries() const override
//    {

//    }

//    virtual void ProcessEntry(const long long entry) override
//    {

//    }


};


//AUTO_REGISTER_PLOTTER(Omega_EpEm_plot)
