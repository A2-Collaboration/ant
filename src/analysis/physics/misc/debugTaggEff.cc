#include "debugTaggEff.h"

#include "slowcontrol/SlowControlVariables.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

using namespace std;


size_t debugTaggEff::getNchannels()
{
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");
    return Tagger->GetNChannels();
}

debugTaggEff::debugTaggEff(const string& name, OptionsPtr opts):
    Physics(name, opts),
    nchannels(getNchannels())
{
    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::Beam->Request();
    slowcontrol::Variables::TaggEff->Request();

    TaggEffTree.CreateBranches(HistFac.makeTTree("taggEff"));

    TaggEffTree.TaggEffs().resize(nchannels);
    TaggEffTree.TaggEffErrors().resize(nchannels);
}

debugTaggEff::~debugTaggEff() {};

void debugTaggEff::ProcessEvent(const TEvent&, manager_t&)
{
    SeenEvents++;

    if(slowcontrol::Variables::TaggerScalers->HasChanged())
    {
        for (auto i = 0u ; i < nchannels ; ++i)
        {
            const auto taggeff = slowcontrol::Variables::TaggEff->Get(i);
            LOG(INFO) << taggeff.Value << "    " << taggeff.Error;
            TaggEffTree.TaggEffs().at(i) = taggeff.Value;
            TaggEffTree.TaggEffErrors().at(i) = taggeff.Error;
        }
        TaggEffTree.TaggerOrRate = slowcontrol::Variables::TaggerScalers->GetTaggerOr();
        TaggEffTree.P2Rate       = slowcontrol::Variables::Beam->GetIonChamber();
        LOG(INFO) << " ladder/p2 = "
                  << TaggEffTree.TaggerOrRate << " / " << TaggEffTree.P2Rate
                  << " = " << 1.0 * TaggEffTree.TaggerOrRate / TaggEffTree.P2Rate;
        TaggEffTree.Tree->Fill();
    }
}

void debugTaggEff::Finish()
{
}

void debugTaggEff::ShowResult()
{
    canvas c("TaggEffs");
    for (auto i = 0u ; i < nchannels ; ++i )
    {
        c << TTree_drawable(TaggEffTree.Tree,std_ext::formatter() << "TaggEffs[" << i << "]");
    }
    c << endc;

    canvas d("Ladder/P2");
    d << TTree_drawable(TaggEffTree.Tree,"TaggerOrRate")
      << TTree_drawable(TaggEffTree.Tree,"P2Rate")
      << TTree_drawable(TaggEffTree.Tree,"1.0 * TaggerOrRate / P2Rate")
      << endc;
}

AUTO_REGISTER_PHYSICS(debugTaggEff)
