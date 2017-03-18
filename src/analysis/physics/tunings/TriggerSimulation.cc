#include "TriggerSimulation.h"

#include "expconfig/ExpConfig.h"
#include "utils/uncertainties/Interpolated.h"

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
    if(fitters.size()<nPhotons+1) {
        fitters.resize(nPhotons+1);
        fitters[nPhotons] = std_ext::make_unique<utils::KinFitter>(
                    std_ext::formatter() << "Fitter_" << nPhotons,
                    nPhotons, fit_model, true
                    );
    }
    return *fitters[nPhotons];
}

struct ProtonPhotonCombs {
    struct comb_t {
        TParticleList Photons;
        TParticlePtr  Proton;
        // set when filtered
        TTaggerHit TaggerHit{0, std_ext::NaN, std_ext::NaN};
        double MissingMass = std_ext::NaN;
    };
    using Combinations_t = std::list<comb_t>;
    const Combinations_t Combinations;

    ProtonPhotonCombs(const TCandidateList& cands) :
        Combinations(
            // use anonymous lambda to create const combinations from candidates,
            // forces filters to copy the combinations (it's cheap because all are shared_ptr)
            [] (const TCandidateList& cands) {
        Combinations_t combs;
        for(auto cand_proton : cands.get_iter()) {
            combs.emplace_back();
            auto& comb = combs.back();
            comb.Proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);
            for(auto cand_photon : cands.get_iter()) {
                if(cand_photon == cand_proton)
                    continue;
                comb.Photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_photon));
            }
        }
        return combs;
    }(cands))
    {}

    static constexpr IntervalD no_cut{-std_ext::inf, std_ext::inf};

//    Combinations_t Filter(const TTaggerHit& taggerhit, const IntervalD& missingmass_cut = {std_ext::, )

};

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

    if(recon.Candidates.size()<3)
        return;

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

        t.nPhotons = recon.Candidates.size()-1; // is at least 2

        // find the fitter for that
//        auto& fitter = GetFitter(t.nPhotons-1);

        t.FitProb = std_ext::NaN;
        // loop over the protons

//            const auto& result = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

//        }


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