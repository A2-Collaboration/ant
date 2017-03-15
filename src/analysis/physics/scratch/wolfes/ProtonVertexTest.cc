#include "ProtonVertexTest.h"

#include "base/vec/LorentzVec.h"
#include "base/Logger.h"

#include "utils/uncertainties/Interpolated.h"

#include "expconfig/setups/Setup.h"

#include "analysis/physics/production/triplePi0.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


ProtonVertexTest::ProtonVertexTest(const string& name, OptionsPtr opts):
     Physics(name, opts),
     uncertModel(utils::UncertaintyModels::Interpolated::makeAndLoad()),
     kinFitterEMB("fitterEMB", 6, uncertModel, true)
{
    kinFitterEMB.SetZVertexSigma(fitter_ZVertex);

    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

    promptrandom.AddPromptRange(Range_Prompt);
    for ( const auto& range: Ranges_Random)
        promptrandom.AddRandomRange(range);

    hist_steps = HistFac.makeTH1D("steps","","# evts.",BinSettings(1,0,0),"steps");
    hist_theta = HistFac.makeTH1D("theta","#theta_{p} [#circ]","#",BinSettings(180));

    tree.CreateBranches(HistFac.makeTTree("ptree"));
}

void ProtonVertexTest::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data   = event.Reconstructed();

    FillStep("seen");

    if (cutOn("ncands", NCands, data.Candidates.size())) return;

    const auto cbAvgTime = data.Trigger.CBTiming;

    for ( const auto& taggerHit: data.TaggerHits )
    {
        FillStep("seen taggerhits");


        promptrandom.SetTaggerHit(taggerHit.Time - cbAvgTime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        tree.prW() = promptrandom.FillWeight();
        FillStep("Promt");



    }





}

void ProtonVertexTest::ShowResult()
{
    canvas("summary")
            << hist_steps
            << hist_theta
            << endc;
}
