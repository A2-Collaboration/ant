#include "MCChannels.h"
#include "base/ParticleTypeTree.h"
#include "analysis/utils/ParticleTools.h"
#include "analysis/plot/RootDraw.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCChannels::MCChannels(const string &name, OptionsPtr opts):
    Physics(name, opts)
{

    const auto prepare_channel_axis = [] (TAxis* axis) {
        for(size_t i=0;i<ParticleTypeTreeDatabase::NumChannels(); ++i) {
            const auto ch = static_cast<ParticleTypeTreeDatabase::Channel>(i);
            const auto str = utils::ParticleTools::GetDecayString(ParticleTypeTreeDatabase::Get(ch), true);
            axis->SetBinLabel(i+1, str.c_str());
        }
        axis->SetBinLabel(axis->GetNbins()-2, "Unknown");
        axis->SetBinLabel(axis->GetNbins()-1, "No tree");
        axis->SetBinLabel(axis->GetNbins(),   "Sum");
    };

    const AxisSettings axis_numChannels("", BinSettings(ParticleTypeTreeDatabase::NumChannels()+3));

    h_database = HistFac.makeTH1D("Channels (in database)", axis_numChannels,
                                  "h_database");
    prepare_channel_axis(h_database->GetXaxis());

    try {
        auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
        const AxisSettings axis_TaggCh(string(Detector_t::ToString(Tagger->Type))+" Channel", {Tagger->GetNChannels()});
        h_database_taggch = HistFac.makeTH2D("Channels (in database) per tagger channel",
                                             axis_numChannels, axis_TaggCh,
                                             "h_database_taggch");
        prepare_channel_axis(h_database_taggch->GetXaxis());
    }
    catch(ExpConfig::ExceptionNoSetup&) {
        LOG(WARNING) << "Disabled per tagger channel histogram as no tagger found in setup";
    }
}

MCChannels::~MCChannels()
{}

void MCChannels::ProcessEvent(const TEvent& event, manager_t&)
{
    total++;

    // Fill is offseted by 1 w.r.t. to bin numbering
    h_database->Fill(h_database->GetNbinsX()-1);
    if(h_database_taggch) {
        for(auto& taggerhit : event.MCTrue().TaggerHits)
            h_database_taggch->Fill(h_database_taggch->GetNbinsX()-1, taggerhit.Channel);
    }
    const auto& ptree = event.MCTrue().ParticleTree;

    if(ptree) {
        counter_production[utils::ParticleTools::GetProductionChannelString(ptree)]++;
        ParticleTypeTreeDatabase::Channel channel;
        if(utils::ParticleTools::TryFindParticleDatabaseChannel(ptree,channel)) {
            h_database->Fill(static_cast<int>(channel));
            if(h_database_taggch) {
                for(auto& taggerhit : event.MCTrue().TaggerHits)
                    h_database_taggch->Fill(static_cast<int>(channel), taggerhit.Channel);
            }
        }
        else {
            h_database->Fill(h_database->GetNbinsX()-3);
            if(h_database_taggch) {
                for(auto& taggerhit : event.MCTrue().TaggerHits)
                    h_database_taggch->Fill(h_database_taggch->GetNbinsX()-3, taggerhit.Channel);
            }
        }
    } else {
        noTree++;
        // Fill is offseted by 1 w.r.t. to bin numbering
        h_database->Fill(h_database->GetNbinsX()-2);
        if(h_database_taggch) {
            for(auto& taggerhit : event.MCTrue().TaggerHits)
                h_database_taggch->Fill(h_database_taggch->GetNbinsX()-2, taggerhit.Channel);
        }
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
    canvas c(GetName());
    c << h_production << h_database;
    if(h_database_taggch)
        c << drawoption("colz") << h_database_taggch;
    c << endc;
}

AUTO_REGISTER_PHYSICS(MCChannels)

