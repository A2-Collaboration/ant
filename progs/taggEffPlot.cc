#include "detail/taggEffClasses.h"

#include "tclap/CmdLine.h"
#include "base/Logger.h"
#include "base/std_ext/system.h"
#include "base/std_ext/string.h"

#include "analysis/plot/RootDraw.h"

#include "TGraph.h"
#include "TRint.h"
#include "TF1.h"

#include <algorithm>

using namespace ant;
using namespace std;
using namespace ant::progs::taggeff;


static TGraph* graph2d(nullptr);
void initGraph(const string& xtitle, const string& ytitle)
{
    graph2d = new TGraph();
    graph2d->GetXaxis()->SetTitle(xtitle.c_str());
    graph2d->GetYaxis()->SetTitle(ytitle.c_str());
    string title = std_ext::formatter() << xtitle << " vs " << ytitle;
    graph2d->SetTitle(title.c_str());
}


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

    auto plot_Lt       = cmd.add<TCLAP::SwitchArg>("","livetime", "plot livetime for bkg1 and bkg2");
    auto plot_rate     = cmd.add<TCLAP::SwitchArg>("","rate", "plot rate for bkg1 and bkg2");
    auto modes = {&plot_Lt,&plot_rate};

    auto cmd_batchmode   = cmd.add<TCLAP::SwitchArg>("b","batch",  "Run in batch mode (no ROOT shell afterwards)");

    auto cmd_filelist   = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","inputfiles to read from",true,"inputfiles");

    cmd.parse(argc,argv);

    auto inputCount = 0u;
    for ( auto m: modes )
        inputCount+=m->get()->isSet();
    if (inputCount!=1)
    {
        string msg = "Exactly one mode is allowed:  ";
        for ( auto m: modes)
            msg += std_ext::formatter() << "--" << m->get()->getName() << "  ";
        failExit(msg);
    }
    auto fileList = cmd_filelist->getValue();

    if (plot_Lt->isSet())
    {
        if (fileList.size() < 1)
            failExit("Provide files");
        plotLiveTimes(fileList);
    }
    if (plot_rate->isSet())
    {
        if (fileList.size() < 1)
            failExit("Provide files");
        plotRates(fileList);
    }

    argc=1; // prevent TRint to parse any cmdline except prog name

    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("plots",&argc,argv,nullptr,0,true);

    if(app && graph2d) {
        auto ytitle = plot_rate->isSet() ? "rate [Hz]" : "livetime [%]";
        graph2d->GetXaxis()->SetTitle("time [s]");
        graph2d->GetYaxis()->SetTitle(ytitle);
        canvas("TaggEff")
         << drawoption("AP")
         << graph2d
         << endc;
        app->Run(kTRUE); // really important to return...
    }


    return EXIT_SUCCESS;
}



void plotLiveTimes(const vector<string>& fileList)
{
    list<treeLoader_t*> tContainers;
    for (const auto& fN: fileList)
        tContainers.emplace_back(new treeLoader_t(fN));

    if ( graph2d )
        throw runtime_error("graph already exists, this should never happen");

    graph2d = timedData::getLtVsTime(tContainers,analysis::HistogramFactory("taggEffPlot"));
}

void plotRates(const vector<string>& fileList)
{
    list<treeLoader_t*> tContainers;
    for (const auto& fN: fileList)
        tContainers.emplace_back(new treeLoader_t(fN));

    if ( graph2d )
        throw runtime_error("graph already exists, this should never happen");

    graph2d = timedData::getRatesVsTime(tContainers,analysis::HistogramFactory("taggEffPlot"));

}
