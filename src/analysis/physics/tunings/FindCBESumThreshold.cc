#include "FindCBESumThreshold.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;




FindCBESumThreshold::FindCBESumThreshold(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(*ExpConfig::Setup::GetLastFound()),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad()),
    fitter("KinFit", 4, fit_model, true) // enable Z vertex by default
{
    fitter.SetZVertexSigma(0); // use unmeasured z vertex

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");
}

void FindCBESumThreshold::ProcessEvent(const TEvent& event, manager_t&)
{
    const TEventData& data = event.Reconstructed();

    steps->Fill("Seen",1);

    const auto& cands = data.Candidates;
    if(cands.size() != 5)
        return;
    steps->Fill("nCands==5",1);

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;


    }
}

void FindCBESumThreshold::ShowResult()
{
    canvas(GetName())
            << steps
            << endc;
}

AUTO_REGISTER_PHYSICS(FindCBESumThreshold)