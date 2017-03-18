#include "TriggerSimulation.h"

#include "utils/uncertainties/Interpolated.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

TriggerSimulation::TriggerSimulation(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(*ExpConfig::Setup::GetLastFound()),
    Clusters_All(HistogramFactory("Clusters_All",HistFac,"Clusters_All")),
    Clusters_Tail(HistogramFactory("Clusters_Tail",HistFac,"Clusters_Tail")),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad()),
    fitter("KinFit", 4, fit_model, true) // enable Z vertex by default
{
    fitter.SetZVertexSigma(0); // use unmeasured z vertex

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

    const AxisSettings axis_CBESum{"CBESum / MeV", {1600, 0, 1600}};
    const AxisSettings axis_CBTiming("CB Timing / ns",{300,-15,10});

    h_CBESum_raw = HistFac.makeTH1D("CBESum raw ",axis_CBESum,"h_CBESum_raw");
    h_CBESum_pr  = HistFac.makeTH1D("CBESum raw prompt-random subtracted",axis_CBESum,"h_CBESum_pr");

    h_CBTiming       = HistFac.makeTH1D("CB Timing", axis_CBTiming, "h_CBTiming");
    h_CBTiming_CaloE = HistFac.makeTH2D("CB Timing vs. CaloE",axis_CBTiming,{"CaloE / MeV", {200,0,100}},"h_CBTiming_CaloE");

    const BinSettings bins_TaggT(200,-30,30);
    h_TaggT = HistFac.makeTH1D("Tagger Timing",{"t_{Tagger}", bins_TaggT},"h_TaggT");
    h_TaggT_corr = HistFac.makeTH1D("Tagger Timing Corrected",{"t_{Tagger} Corrected", bins_TaggT},"h_TaggT_corr");
    h_TaggT_CBTiming = HistFac.makeTH2D("Tagger Timing vs. CBTiming",{"t_{Tagger}", bins_TaggT},axis_CBTiming,"h_TaggT_CBTiming");
}

TriggerSimulation::ClusterPlots_t::ClusterPlots_t(const HistogramFactory& HistFac)
{
    const AxisSettings axis_CaloE("CaloE / MeV",{100,0,20});
    const AxisSettings axis_ClSize("ClusterSize",{10});
    const AxisSettings axis_nCl("nClusters",{10});
    const AxisSettings axis_timing("t / ns",{100,-30,30});

    h_CaloE_ClSize = HistFac.makeTH2D("CaloE vs. ClusterSize", axis_CaloE, axis_ClSize, "h_CaloE_ClSize");
    h_CaloE_nCl = HistFac.makeTH2D("CaloE vs. nClusters", axis_CaloE, axis_nCl, "h_CaloE_nCl");
    h_CaloE_Time = HistFac.makeTH2D("CaloE vs. Time",axis_CaloE, axis_timing,"h_CaloE_Time");
    h_Hits_stat = HistFac.makeTH1D("Hits status","","",BinSettings(4),"h_Hits_stat");
    h_Hits_E_t  = HistFac.makeTH2D("ClHits Energy vs. Time",{"E_{hit} / MeV",{100,0,50}}, axis_timing ,"h_Hits_E_t");
}

void TriggerSimulation::ClusterPlots_t::Fill(const TEventData& recon) const
{
    for(const TCluster& cluster : recon.Clusters) {
        if(cluster.DetectorType == Detector_t::Type_t::CB) {
            h_CaloE_ClSize->Fill(cluster.Energy,cluster.Hits.size());
            h_CaloE_nCl->Fill(cluster.Energy,recon.Clusters.size());
            h_CaloE_Time->Fill(cluster.Energy, cluster.Time);
            for(const auto& hit : cluster.Hits) {
                h_Hits_E_t->Fill(hit.Energy, hit.Time);
                h_Hits_stat->Fill("Seen",1.0);
                if(hit.IsSane())
                    h_Hits_stat->Fill("Sane",1.0);
                if(isfinite(hit.Time))
                    h_Hits_stat->Fill("Time ok",1.0);
                if(isfinite(hit.Time))
                    h_Hits_stat->Fill("Energy ok",1.0);
            }
        }
    }
}

void TriggerSimulation::ClusterPlots_t::Show(canvas &c) const
{
    c << drawoption("colz")
      << h_CaloE_ClSize << h_CaloE_nCl << h_CaloE_Time
      << h_Hits_stat << h_Hits_E_t
      << endr;
}

void TriggerSimulation::ProcessEvent(const TEvent& event, manager_t&)
{
    steps->Fill("Seen",1);

    if(!triggersimu.ProcessEvent(event)) {
        steps->Fill("TriggerSimu failed", 1.0);
        return;
    }

    steps->Fill("Triggered", triggersimu.HasTriggered());

    const auto& recon = event.Reconstructed();

    h_CBESum_raw->Fill(triggersimu.GetCBEnergySum());
    h_CBTiming->Fill(triggersimu.GetRefTiming());
    for(const TCluster& cluster : recon.Clusters) {
        if(cluster.DetectorType == Detector_t::Type_t::CB) {
            h_CBTiming_CaloE->Fill(triggersimu.GetRefTiming(),cluster.Energy);
        }
    }

    Clusters_All.Fill(recon);
    if(IntervalD(-10,-5).Contains(triggersimu.GetRefTiming())) {
        // investigate the tail
        Clusters_Tail.Fill(recon);
    }


    for(const TTaggerHit& taggerhit : recon.TaggerHits) {

        steps->Fill("Seen taggerhits",1.0);

        h_TaggT->Fill(taggerhit.Time);
        h_TaggT_CBTiming->Fill(taggerhit.Time, triggersimu.GetRefTiming());
        const auto& taggertime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        h_TaggT_corr->Fill(taggertime);

        promptrandom.SetTaggerTime(taggertime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        steps->Fill("Acc taggerhits",1.0);

        h_CBESum_pr->Fill(triggersimu.GetCBEnergySum(), promptrandom.FillWeight());

        const auto& cands = recon.Candidates;
        if(cands.size() != 5)
            return;
        steps->Fill("nCands==5",1);
    }
}

void TriggerSimulation::ShowResult()
{
    canvas(GetName()) << drawoption("colz")
            << steps
            << h_TaggT << h_TaggT_CBTiming << h_TaggT_corr
            << h_CBTiming
            << h_CBESum_raw << h_CBESum_pr
            << endc;
    canvas c(GetName()+": CBTiming Tail");
    Clusters_All.Show(c);
    Clusters_Tail.Show(c);
    c << endc;
}





AUTO_REGISTER_PHYSICS(TriggerSimulation)