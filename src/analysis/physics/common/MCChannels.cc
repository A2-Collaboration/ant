#include "MCChannels.h"
#include "base/ParticleTypeTree.h"
#include "analysis/utils/ParticleTools.h"
#include "analysis/plot/RootDraw.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCChannels::MCChannels(const string &name, OptionsPtr opts):
    Physics(name, opts)
{
    h_channels = HistFac.makeTH1D("Channels (database)", "Channel", "#", BinSettings(ParticleTypeTreeDatabase::NumChannels()+2), "all_channels");
    for(size_t i=0;i<ParticleTypeTreeDatabase::NumChannels(); ++i) {
        const auto ch = static_cast<ParticleTypeTreeDatabase::Channel>(i);
        const auto str = utils::ParticleTools::GetDecayString(ParticleTypeTreeDatabase::Get(ch), true);
        h_channels->GetXaxis()->SetBinLabel(i+1, str.c_str());
    }
    h_channels->GetXaxis()->SetBinLabel(h_channels->GetNbinsX()-1, "No Tree into");
    h_channels->GetXaxis()->SetBinLabel(h_channels->GetNbinsX(),    "Sum");
}

MCChannels::~MCChannels()
{}

void MCChannels::ProcessEvent(const TEvent &event, manager_t&)
{
    total++;
    h_channels->Fill(h_channels->GetNbinsX()-1);
    const auto& head = event.MCTrue().ParticleTree;
    if(head) {
        const auto str = utils::ParticleTools::GetProductionChannelString(head);
        counter[str]++;
        ParticleTypeTreeDatabase::Channel channel;
        ParticleTypeTree typetree;
        if(utils::ParticleTools::TryFindParticleTypeTree(head,channel,typetree)) {
            h_channels->Fill(static_cast<int>(channel));
        }

    } else {
        noTree++;
        h_channels->Fill(h_channels->GetNbinsX()-2);
    }
}
void MCChannels::Finish() {

    hist       = HistFac.makeTH1D("Production Channels", "Channel", "#", BinSettings(2+counter.size()),"channels");

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


}
void MCChannels::ShowResult()
{
    canvas(GetName()) << hist << h_channels << endc;
}

AUTO_REGISTER_PHYSICS(MCChannels)

