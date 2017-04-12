#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"

#include "base/Logger.h"

#include "base/std_ext/system.h"
#include "base/std_ext/string.h"

#include "base/PlotExt.h"
#include "analysis/plot/RootDraw.h"

#include "expconfig/ExpConfig.h"

#include "TGraph.h"
#include "TAxis.h"
#include "TRint.h"

using namespace ant;
using namespace std;


auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};

void plotLiveTimes(const vector<string>& fileList);
void plotRates(const vector<string>& fileList);


int main( int argc, char** argv )

{
    SetupLogger();
    TCLAP::CmdLine cmd("taggEffPlot: plot livetime or rate vs time in seconds from first ProcessTaggEff-output of given list", ' ', "0.1");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",true,"", &allowedsetupnames);

    cmd.parse(argc,argv);

    auto graph = new TGraph();



    const auto& setupname = cmd_setup->getValue();
    try {
        ExpConfig::Setup::SetByName(setupname);
        auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
        for(unsigned ch=0;ch<tagger->GetNChannels();ch++) {
            GraphExt::FillGraph(graph,ch,tagger->GetPhotonEnergy(ch));
        }
    }
    catch(ExpConfig::Exception e) {
        failExit(std_ext::formatter() << "Specified setup '" << setupname << "' did not provide a tagger: " << e.what());
    }

    graph->GetXaxis()->SetTitle("channel");
    graph->GetYaxis()->SetTitle("E_{#gamma} [MeV]");
    graph->SetMarkerStyle(7);

    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = std_ext::make_unique<TRint>("Tagger",&argc,argv,nullptr,0,true);
    canvas("Tagger") << drawoption("ap") << graph << endc;
    app->Run(kTRUE);
    return EXIT_SUCCESS;
}



