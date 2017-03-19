#include "TriggerSimulation.h"

#include "expconfig/ExpConfig.h"
#include "utils/uncertainties/Interpolated.h"
#include "utils/ProtonPhotonCombs.h"
#include "utils/Combinatorics.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

TriggerSimulation::TriggerSimulation(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()),
    Clusters_All(HistogramFactory("Clusters_All",HistFac,"Clusters_All")),
    Clusters_Tail(HistogramFactory("Clusters_Tail",HistFac,"Clusters_Tail")),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad())
{

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(10),"steps");

    const AxisSettings axis_CBESum{"CBESum / MeV", {1600, 0, 1600}};
    const AxisSettings axis_CBTiming("CB Timing / ns",{300,-15,10});

    h_CBESum_raw = HistFac.makeTH1D("CBESum raw ",axis_CBESum,"h_CBESum_raw");
    h_CBESum_pr  = HistFac.makeTH1D("CBESum raw p-r sub",axis_CBESum,"h_CBESum_pr");
    h_CBESum_fit  = HistFac.makeTH1D("CBESum fit p-r sub",axis_CBESum,"h_CBESum_fit");

    h_CBTiming       = HistFac.makeTH1D("CB Timing", axis_CBTiming, "h_CBTiming");
    h_CBTiming_CaloE = HistFac.makeTH2D("CB Timing vs. CaloE",axis_CBTiming,{"CaloE / MeV", {200,0,100}},"h_CBTiming_CaloE");

    const BinSettings bins_TaggT(200,-30,30);
    h_TaggT = HistFac.makeTH1D("Tagger Timing",{"t_{Tagger}", bins_TaggT},"h_TaggT");
    h_TaggT_corr = HistFac.makeTH1D("Tagger Timing Corrected",{"t_{Tagger} Corrected", bins_TaggT},"h_TaggT_corr");
    h_TaggT_CBTiming = HistFac.makeTH2D("Tagger Timing vs. CBTiming",{"t_{Tagger}", bins_TaggT},axis_CBTiming,"h_TaggT_CBTiming");

    t.CreateBranches(HistFac.makeTTree("tree"));
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

utils::KinFitter& TriggerSimulation::GetFitter(unsigned nPhotons)
{
    // lazy init the fitters on demand
    if(fitters.size()<nPhotons+1 || !fitters[nPhotons]) {
        fitters.resize(nPhotons+1);
        fitters[nPhotons] = std_ext::make_unique<utils::KinFitter>(
                    std_ext::formatter() << "Fitter_" << nPhotons,
                    nPhotons, fit_model, true
                    );
        fitters[nPhotons]->SetZVertexSigma(0); // unmeasured z vertex
    }
    return *fitters[nPhotons];
}



void TriggerSimulation::ProcessEvent(const TEvent& event, manager_t&)
{

    steps->Fill("Seen",1);

    if(!triggersimu.ProcessEvent(event)) {
        steps->Fill("TriggerSimu failed", 1.0);
        return;
    }

    steps->Fill("Triggered", triggersimu.HasTriggered());

    // as MC may have also some pure TAPS events,
    // zero CBEnergySum can be suppressed,
    // as we can't return when we're not triggered
    // (that's what we want to determine)
    if(triggersimu.GetCBEnergySum()==0)
        return;

    const auto& recon = event.Reconstructed();

    // gather tree stuff already here, before we forget it :)
    t.IsMC = recon.ID.isSet(TID::Flags_t::MC);
    t.Triggered = triggersimu.HasTriggered();
    t.CBEnergySum = triggersimu.GetCBEnergySum();

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

    // want at least proton and two gammas in final state
    // beyond this point
    if(recon.Candidates.size()<3)
        return;

    utils::ProtonPhotonCombs proton_photons(recon.Candidates);

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

        t.TaggW = promptrandom.FillWeight();
        t.TaggT = taggertime;
        t.TaggE = taggerhit.PhotonEnergy;
        t.TaggCh = taggerhit.Channel;

        t.nPhotons = recon.Candidates.size()-1; // is at least 2

        // setup a very inclusive filter, just to speed up fitting
        auto filtered_proton_photons = proton_photons([this] (const string& cut) { steps->Fill(cut.c_str(), 1.0); }).
                                       FilterIM(). // no filter for IM of photons
                                       FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(500).Round());

        if(filtered_proton_photons.empty()) {
            steps->Fill("No combs left",1.0);
            continue;
        }

        auto& fitter = GetFitter(t.nPhotons);

        // loop over the (filtered) proton combinations
        t.FitProb = std_ext::NaN;
        for(const auto& comb : filtered_proton_photons) {

            const auto& result = fitter.DoFit(comb.TaggerHit.PhotonEnergy, comb.Proton, comb.Photons);

            if(result.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(t.FitProb, result.Probability))
                continue;

            // do combinatorics
            const auto fill_IM_Combs = [] (vector<double>& v, const TParticleList& photons) {
                auto combs = utils::makeCombination(photons, 2);
                v.resize(combs.size());
                for(auto& im : v) {
                    im = (*combs.at(0) + *combs.at(1)).M();
                    combs.next();
                }
            };

            fill_IM_Combs(t.IM_Combs_fitted, fitter.GetFittedPhotons());
            fill_IM_Combs(t.IM_Combs_raw, comb.Photons);
        }

        if(t.FitProb>0.01) {
            steps->Fill("FitProb>0.01",1.0);
            t.Tree->Fill();
            h_CBESum_fit->Fill(t.CBEnergySum, t.TaggW);
        }

    }
}

void TriggerSimulation::ShowResult()
{
    canvas(GetName()) << drawoption("colz")
            << steps
            << h_TaggT << h_TaggT_CBTiming << h_TaggT_corr
            << h_CBTiming
            << h_CBESum_raw << h_CBESum_pr << h_CBESum_fit
            << endc;
    canvas c(GetName()+": CBTiming Tail");
    Clusters_All.Show(c);
    Clusters_Tail.Show(c);
    c << endc;
}





AUTO_REGISTER_PHYSICS(TriggerSimulation)