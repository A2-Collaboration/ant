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
    h_database = HistFac.makeTH1D("Channels (in database)", "", "",
                                         BinSettings(ParticleTypeTreeDatabase::NumChannels()+2), "h_database");
    for(size_t i=0;i<ParticleTypeTreeDatabase::NumChannels(); ++i) {
        const auto ch = static_cast<ParticleTypeTreeDatabase::Channel>(i);
        const auto str = utils::ParticleTools::GetDecayString(ParticleTypeTreeDatabase::Get(ch), true);
        h_database->GetXaxis()->SetBinLabel(i+1, str.c_str());
    }
    h_database->GetXaxis()->SetBinLabel(h_database->GetNbinsX()-1, "No Tree into");
    h_database->GetXaxis()->SetBinLabel(h_database->GetNbinsX(),   "Sum");
}

MCChannels::~MCChannels()
{}

void MCChannels::ProcessEvent(const TEvent &event, manager_t&)
{
    total++;
    h_database->Fill(h_database->GetNbinsX()-1); // Fill is offseted by 1 w.r.t. to bin numbering
    const auto& ptree = event.MCTrue().ParticleTree;
    if(ptree) {
        counter_production[utils::ParticleTools::GetProductionChannelString(ptree)]++;
        ParticleTypeTreeDatabase::Channel channel;
        if(utils::ParticleTools::TryFindParticleDatabaseChannel(ptree,channel)) {
            h_database->Fill(static_cast<int>(channel));
        }
    } else {
        noTree++;
        h_database->Fill(h_database->GetNbinsX()-2); // Fill is offseted by 1 w.r.t. to bin numbering
    }
}
void MCChannels::Finish() {

    h_production = HistFac.makeTH1D("Production Channels", "", "",
                                           BinSettings(2+counter_production.size()),"h_production");

    h_production->SetBinContent(1, total);
    h_production->GetXaxis()->SetBinLabel(1, "Total");

    h_production->SetBinContent(2, noTree);
    h_production->GetXaxis()->SetBinLabel(2, "no tree");

    int b=3;
    for(const auto entry : counter_production) {
        h_production->SetBinContent(b, entry.second);
        h_production->GetXaxis()->SetBinLabel(b, entry.first.c_str());
        ++b;
    }


}
void MCChannels::ShowResult()
{
    canvas(GetName()) << h_production << h_database << endc;
}

AUTO_REGISTER_PHYSICS(MCChannels)

