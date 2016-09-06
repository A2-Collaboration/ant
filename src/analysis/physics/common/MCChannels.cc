#include "MCChannels.h"
#include "base/ParticleTypeTree.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCChannels::MCChannels(const string &name, OptionsPtr opts):
    Physics(name, opts)
{}

MCChannels::~MCChannels()
{}

void MCChannels::ProcessEvent(const TEvent &event, manager_t&)
{
    total++;
    if(event.HasMCTrue()) {
        const auto& head = event.MCTrue().ParticleTree;
        if(head) {
            const auto str = utils::ParticleTools::GetProductionChannelString(head);
            counter[str]++;
        } else {
            noTree++;
        }
    } else {
        noTree++;
    }
}

void MCChannels::ShowResult()
{
    hist=HistFac.makeTH1D("Production Channels", "Channel", "#", BinSettings(2+counter.size()),"channels");

    hist->SetBinContent(1, total);
    hist->GetXaxis()->SetBinLabel(1, "Total");

    hist->SetBinContent(2, noTree);
    hist->GetXaxis()->SetBinLabel(2, "no tree");

    int b=3;
    for(const auto entry : counter) {
        hist->SetBinContent(b, entry.second);
        hist->GetXaxis()->SetBinLabel(b, entry.first.c_str());
        ++b;
    }

    canvas(GetName()) << hist << endc;
}

AUTO_REGISTER_PHYSICS(MCChannels)

