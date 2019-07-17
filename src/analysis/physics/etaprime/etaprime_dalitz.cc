#include "etaprime_dalitz.h"

#include "utils/Combinatorics.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"
#include "analysis/utils/uncertainties/MeasuredProton.h"
#include "base/std_ext/vector.h"
#include "base/std_ext/string.h"
#include "base/std_ext/container.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/TAPS.h"

#include "base/Logger.h"


using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

// make the linker happy
std::once_flag EtapDalitz::Settings_t::initialized;


void EtapDalitz::set_beamtime(common_tree *t)
{
    if (std_ext::contains(ExpConfig::Setup::Get().GetName(), "2014_07"))
        t->beamtime = 1;
    else if (std_ext::contains(ExpConfig::Setup::Get().GetName(), "2014_10"))
        t->beamtime = 2;
    else if (std_ext::contains(ExpConfig::Setup::Get().GetName(), "2014_12"))
        t->beamtime = 3;
}

double EtapDalitzTools::effective_radius(const TCandidatePtr cand) const
{
    return clustertools.EffectiveRadius(*(cand->FindCaloCluster()));
}

double EtapDalitzTools::lat_moment(const TCandidatePtr cand) const
{
    return clustertools.LateralMoment(*(cand->FindCaloCluster()));
}

ParticleTypeTree EtapDalitzTools::base_tree()
{
    ParticleTypeTree t = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(ParticleTypeDatabase::BeamProton);
    t->CreateDaughter(ParticleTypeDatabase::Proton);
    return t;
}

ParticleTypeTree EtapDalitzTools::etap_3g()
{
    auto t = base_tree();
    auto etap = t->CreateDaughter(ParticleTypeDatabase::EtaPrime);
    etap->CreateDaughter(ParticleTypeDatabase::Photon);
    etap->CreateDaughter(ParticleTypeDatabase::Photon);
    etap->CreateDaughter(ParticleTypeDatabase::Photon);
    return t;
}

void EtapDalitzTools::fake_comb_t::reset()
{
    Proton = nullptr;
    Photons.resize(0);
    PhotonSum = LorentzVec();
    MissingMass = std_ext::NaN;
    DiscardedEk = 0.;
}

void EtapDalitzTools::fake_comb_t::calc_values(const TTaggerHit& taggerhit)
{
    for (const auto& p : Photons)
        PhotonSum += *p;

    const auto beam_target = taggerhit.GetPhotonBeam() + LorentzVec::AtRest(ParticleTypeDatabase::Proton.Mass());
    MissingMass = (beam_target - PhotonSum).M();
}

APLCON::Fit_Settings_t EtapDalitz::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    return settings;
}

EtapDalitz::PerChannel_t::PerChannel_t(const std::string& Name, const string& Title, HistogramFactory& hf):
    title(Title),
    name(Name)
{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);
    const BinSettings pid_channels(detector->GetNChannels());
    const BinSettings veto_energy(1000, 0, 10);
    const BinSettings energy(1200);

    steps = hf.makeTH1D(title + " Accepted Events", "step", "#", BinSettings(10), name + " steps");

    etapIM = hf.makeTH1D(title + " IM #eta' all comb", "IM [MeV]", "#", energy, name + " etapIM");
    MM = hf.makeTH1D(title + " Missing Mass proton", "MM [MeV]", "#", BinSettings(1600), name + " MM");
    trueZVertex = hf.makeTH1D(title + " true Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " trueZ");

    treefitProb = hf.makeTH1D(title + " treefitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " treefitProb");
    treefit_freeZ_prob = hf.makeTH1D(title + " free Z treefitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " treefit_freeZ_prob");
    kinfitProb = hf.makeTH1D(title + " kinfitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " kinfitProb");
    kinfit_freeZ_prob = hf.makeTH1D(title + " free Z kinfitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " kinfit_freeZ_prob");
    kinfit_ZVertex = hf.makeTH1D(title + " kinfitted Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " kinfit_ZVertex");
    kinfit_freeZ_ZVertex = hf.makeTH1D(title + " free Z kinfitted Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " kinfit_freeZ_ZVertex");
    treefit_ZVertex = hf.makeTH1D(title + " treefitted Z Vertex", "z [cm]", "#", BinSettings(300, -30, 30), name + " treefit_ZVertex");
    treefit_freeZ_ZVertex = hf.makeTH1D(title + " free Z treefitted Z Vertex", "z [cm]", "#", BinSettings(300, -30, 30), name + " treefit_freeZ_ZVertex");
    antiPionProb = hf.makeTH1D(title + " Probability anti-#pi Fit", "probability", "#", BinSettings(500, 0, 1), name + " antiPionProb");

    if (Settings_t::get().less_plots())
        return;

    etapIM_final = hf.makeTH1D(title + " IM #eta' final", "IM [MeV]", "#", energy, name + " etapIM_final");
    IM2d = hf.makeTH2D(title + " IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), name + " IM2d");

    effect_rad = hf.makeTH1D(title + " Effective Radius", "R", "#", BinSettings(500, 0, 50), name + " effect_rad");
    effect_rad_E = hf.makeTH2D(title + " Effective Radius vs. Cluster Energy", "E [MeV]", "R", energy, BinSettings(500, 0, 50), name + " effect_rad_E");
    cluster_size = hf.makeTH1D(title + " Cluster Size", "N", "#", BinSettings(50), name + " cluster_size");
    cluster_size_E = hf.makeTH2D(title + " Cluster Size vs. Cluster Energy", "E [MeV]", "N", energy, BinSettings(50), name + " cluster_size_E");
    lat_moment = hf.makeTH1D(title + " Lateral Moment", "L", "#", BinSettings(200, 0, 1), name + " lat_moment");
    lat_moment_E = hf.makeTH2D(title + " Lateral Moment vs. Cluster Energy", "E [MeV]", "L", energy, BinSettings(200, 0, 1), name + " lat_moment_E");

    proton_E_theta = hf.makeTH2D(title + " proton", "E [MeV]", "#vartheta [#circ]", energy, BinSettings(360, 0, 180), name + " e_theta");
}

void EtapDalitz::PerChannel_t::Show()
{
    if (Settings_t::get().less_plots()) {
        canvas("Per Channel: " + title) << steps
                                        << kinfitProb
                                        << endc;
    } else {
        canvas("Per Channel: " + title) << steps
                                        << etapIM
                                        << etapIM_final
                                        << trueZVertex
                                        << kinfitProb
                                        << kinfit_freeZ_prob
                                        << treefitProb
                                        << treefit_freeZ_prob
                                        << endc;
    }
}

void EtapDalitz::PerChannel_t::Fill(const TEventData& d)
{
    if (Settings_t::get().less_plots())
        return;

    auto particles = d.ParticleTree ?
                         utils::ParticleTypeList::Make(d.ParticleTree) :
                         utils::ParticleTypeList::Make(d.Candidates);
    const auto& protons = particles.Get(ParticleTypeDatabase::Proton);
    if (!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), p->Theta()*TMath::RadToDeg());
    }
}

static int getDetectorAsInt(const Detector_t::Any_t& d)
{
    if (d & Detector_t::Type_t::CB)
        return 1;
    else if (d & Detector_t::Type_t::TAPS)
        return 2;

    return 0;
}

EtapDalitz::EtapDalitz(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    model_sergey(make_shared<utils::UncertaintyModels::FitterSergey>()),  // default base model
    model_data(utils::UncertaintyModels::Interpolated::makeAndLoad(
                   utils::UncertaintyModels::Interpolated::Type_t::Data,
                   // use Sergey as starting point
                   model_sergey
                   )),
    model_MC(utils::UncertaintyModels::Interpolated::makeAndLoad(
                 utils::UncertaintyModels::Interpolated::Type_t::MC,
                 // use Sergey as starting point
                 model_sergey
                 )),
    model_data_protonMeasured(utils::UncertaintyModels::Interpolated::makeAndLoad(
                                  utils::UncertaintyModels::Interpolated::Type_t::Data,
                                  // use measured proton uncertainty with Sergey as base
                                  make_shared<utils::UncertaintyModels::MeasuredProton>(model_sergey),
                                  true)),  // use proton energy
    model_MC_protonMeasured(utils::UncertaintyModels::Interpolated::makeAndLoad(
                                utils::UncertaintyModels::Interpolated::Type_t::MC,
                                // use measured proton uncertainty with Sergey as base
                                make_shared<utils::UncertaintyModels::MeasuredProton>(model_sergey),
                                true)),  // use proton energy
    kinfit(nullptr, opts->HasOption("SigmaZ"), MakeFitSettings(20)),
    kinfit_freeZ(nullptr, true,                MakeFitSettings(20)),
    treefitter_etap(etap_3g(), nullptr,
                    opts->HasOption("SigmaZ"), {}, MakeFitSettings(20)
                    ),
    treefitter_etap_freeZ(etap_3g(), nullptr,
                          true, {}, MakeFitSettings(20)
                          )
{
    //promptrandom.AddPromptRange({-5, 5});
    promptrandom.AddPromptRange({-3, 2});
    promptrandom.AddRandomRange({-35, -10});
    promptrandom.AddRandomRange({10, 35});

    // initialize settings
    settings.init(opts->Get<bool>("reference", 0),
                  opts->Get<bool>("reference_only", 0),
                  opts->Get<bool>("less_plots", 0));

    cb = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);
    taps = make_shared<expconfig::detector::TAPS_2013_11>(false, false, false);

    sig.CreateBranches(HistFac.makeTTree("signal"));
    if (settings.reference()) {
        if (settings.reference_only())
            LOG(INFO) << "Only Reference channel will be analysed";
        else
            LOG(INFO) << "Reference channel included in analysis";
        ref.CreateBranches(HistFac.makeTTree("ref"));
        etap2g = new Etap2g("Etap2g", opts);
        etap2g->linkTree(ref);
        etap2g->setPromptRandom(promptrandom);
    }

    if (settings.less_plots())
        LOG(INFO) << "Less histograms will be created and stored";

    const BinSettings tagger_time_bins(2000, -200, 200);

    h_tagger_time = HistFac.makeTH1D("Tagger Time", "t [ns]", "#", tagger_time_bins, "h_tagger_time");
    h_tagger_time_CBavg = HistFac.makeTH1D("Tagger Time - CB avg time", "t [ns]", "#", tagger_time_bins, "h_tagger_time_CBavg");

    h_counts = HistFac.makeTH1D("Events per Channel", "channel", "#", BinSettings(20), "h_counts");
    h_nCands = HistFac.makeTH1D("Number of Candidates", "#Candidates", "#", BinSettings(30), "h_nCands");
    missed_channels = HistFac.makeTH1D("Unlisted Channels", "", "Total Events seen", BinSettings(20), "missed_channels");
    found_channels  = HistFac.makeTH1D("Listed Channels", "", "Total Events seen", BinSettings(20), "found_channels");

    const BinSettings energybins(1000, 0, 10);

    if (!settings.less_plots()) {
        h_cluster_CB = HistFac.makeTH1D("#Cluster CB", "#Cluster", "#", BinSettings(20), "h_cluster_CB");
        h_cluster_TAPS = HistFac.makeTH1D("#Cluster TAPS", "#Cluster", "#", BinSettings(20), "h_cluster_TAPS");

        h_pTheta = HistFac.makeTH1D("#vartheta proton candidate", "#vartheta_{p} [#circ]", "#", BinSettings(720, 0, 180), "h_pTheta");
        h_protonVeto = HistFac.makeTH1D("Veto energy identified proton", "Veto [MeV]", "#", energybins, "h_protonVeto");
        h_etapIM_final = HistFac.makeTH1D("IM #eta' final", "IM [MeV]", "#", BinSettings(1200), "h_etapIM_final");
        h_IM2d = HistFac.makeTH2D("IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), "h_IM2d");
        h_etap = HistFac.makeTH2D("Kinematics #eta'", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(360, 0, 180), "h_etap");
        h_proton = HistFac.makeTH2D("Kinematics p", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(160, 0, 80), "h_proton");
        h_subIM_2g = HistFac.makeTH1D("#pi^{0} Candidate sub IM 2#gamma", "IM [MeV]", "#", BinSettings(1600, 0, 400), "h_subIM_2g");
        h_subIM_2g_fit = HistFac.makeTH1D("#pi^{0} Candidate sub IM 2#gamma after KinFit", "IM [MeV]", "#", BinSettings(1600, 0, 400), "h_subIM_2g_fit");
    }

    // get target information
    const auto target = ExpConfig::Setup::Get().GetTargetProperties();

    // set sigma to 0 for unmeasured --> free z vertex
    kinfit_freeZ.SetZVertexSigma(0);
    kinfit_freeZ.SetTarget(target.length);
    treefitter_etap_freeZ.SetZVertexSigma(0);
    treefitter_etap_freeZ.SetTarget(target.length);

    if (opts->HasOption("SigmaZ")) {
        double sigma_z = 0.;
        std_ext::copy_if_greater(sigma_z, opts->Get<double>("SigmaZ", 0.));
        LOG(INFO) << "Fit Z vertex enabled with sigma = " << sigma_z;
        kinfit.SetZVertexSigma(sigma_z);
        treefitter_etap.SetZVertexSigma(sigma_z);
    }

    // setup does never change, so set beamtime information once and for all
    set_beamtime(&sig);
    if (settings.reference())
        set_beamtime(&ref);
}

void EtapDalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);
    sig.init();

    const auto& data = event.Reconstructed();
    const bool MC = data.ID.isSet(TID::Flags_t::MC);


    sig.MCtrue = MC;
    sig.channel = reaction_channels.identify(event.MCTrue().ParticleTree);
    if (MC && !sig.channel)  // assign other_index in case of an empty or unknown particle tree for MC (tagged as data otherwise)
        sig.channel = reaction_channels.other_index;
    sig.trueZVertex = event.MCTrue().Target.Vertex.z;  // NaN in case of data

    if (sig.channel == ReactionChannelList_t::other_index) {
        if (MC)
            missed_channels->Fill(utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree).c_str(), 1);
    } else
        found_channels->Fill(sig.channel);

    // identify the currently processed channel
    channel_id(event, chan_id);

    // manage histogram structure for different channels, get histograms for current channel
    auto h = manage_channel_histograms_get_current(MC, event);
    h.trueZVertex->Fill(sig.trueZVertex);

    const auto& cands = data.Candidates;
    //const auto nCandidates = cands.size();
    sig.nCands = cands.size();
    h_nCands->Fill(sig.nCands);
    h.steps->Fill("seen", 1);

    // histogram amount of CB and TAPS clusters
    if (!settings.less_plots()) {
        size_t nCB = 0, nTAPS = 0;
        count_clusters(cands, nCB, nTAPS);
        h_cluster_CB->Fill(nCB);
        h_cluster_TAPS->Fill(nTAPS);
    }

    if (!triggersimu.HasTriggered())
        return;
    h.steps->Fill("triggered", 1);

    sig.CBSumE = triggersimu.GetCBEnergySum();

    sig.CBAvgTime = triggersimu.GetRefTiming();
    if (!isfinite(sig.CBAvgTime))
        return;
    h.steps->Fill("CBAvgTime OK", 1);

    // set fitter uncertainty models
    {
        const auto& model = MC ? model_MC : model_data;
        kinfit.SetUncertaintyModel(model);
        kinfit_freeZ.SetUncertaintyModel(model);
        treefitter_etap.SetUncertaintyModel(model);
        treefitter_etap_freeZ.SetUncertaintyModel(model);
    }

    if (settings.reference()) {
        ref.init();
        ref.MCtrue = sig.MCtrue;
        ref.channel = sig.channel;
        ref.trueZVertex = sig.trueZVertex;
        ref.nCands = sig.nCands;
        ref.CBSumE = sig.CBSumE;
        ref.CBAvgTime = sig.CBAvgTime;

        etap2g->Process(event);

        if (settings.reference_only())
            return;
    }


//    if (cands.size() != Cuts_t::N_FINAL_STATE)
//        return;
//    h.steps->Fill("#cands", 1);

    // q2 preselection on MC data
    if (MC && Cuts_t::Q2_PRESELECTION) {
        if (!q2_preselection(event.MCTrue(), Cuts_t::Q2_MIN_VALUE))
            return;
        stringstream ss;
        ss << "MC q2 > " << Cuts_t::Q2_MIN_VALUE;
        h.steps->Fill(ss.str().c_str(), 1);
    }

    //const auto mass_etap = ParticleTypeDatabase::EtaPrime.Mass();
    //const interval<double> etap_im({mass_etap-Cuts_t::ETAP_SIGMA, mass_etap+Cuts_t::ETAP_SIGMA});

    utils::ProtonPhotonCombs proton_photons(cands);


    double best_prob_fit = -std_ext::inf;
    // loop over all tagger hits
    for (const TTaggerHit& taggerhit : data.TaggerHits) {
        sig.reset();

        if (!MC) {
            h_tagger_time->Fill(taggerhit.Time);
            h_tagger_time_CBavg->Fill(triggersimu.GetCorrectedTaggerTime(taggerhit));
        }

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h.steps->Fill("time window", 1);

        sig.TaggW = promptrandom.FillWeight();
        sig.TaggE = taggerhit.PhotonEnergy;
        sig.TaggT = taggerhit.Time;
        sig.TaggCh = taggerhit.Channel;

//        particle_combs_t selection = proton_photons()
//                .Observe([h] (const std::string& s) { h.steps->Fill(s.c_str(), 1.); }, "[S] ")
//                // require 3 photons and allow discarded energy of 100 MeV
//                .FilterMult(settings.n_final_state_etap, settings.max_discarded_energy)
//                .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(settings.mm_window_size).Round())  // MM cut on proton mass
//                .FilterCustom([=] (const particle_comb_t& p) {
//                    // ensure the possible proton candidate is kinematically allowed
//                    if (std_ext::radian_to_degree(p.Proton->Theta()) > settings.max_proton_theta)
//                        return true;
//                    return false;
//                }, "proton #vartheta")
//                .FilterCustom([] (const particle_comb_t& p) {
//                    // require 2 PID entries for the eta' candidate
//                    if (std::count_if(p.Photons.begin(), p.Photons.end(),
//                                      [](TParticlePtr g){ return g->Candidate->VetoEnergy > .3; }) < 2)
//                        return true;
//                    return false;
//                }, "2 PIDs");

//        if (selection.empty())
//            continue;
//        h.steps->Fill("Selection", 1);

        //begin of test to only use one combination with 3CB and 1TAPS cluster
        size_t nCB = 0, nTAPS = 0;
        count_clusters(cands, nCB, nTAPS);
        if (nCB != 3 || nTAPS != 1)
            return;
        h.steps->Fill("3CB&1TAPS", 1);

        fake_comb_t comb;
        comb.reset();

        for (auto c : cands.get_iter())
            if (c->Detector & Detector_t::Type_t::CB)
                comb.Photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, c));
            else if (c->Detector & Detector_t::Type_t::TAPS)
                comb.Proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, c);

        comb.calc_values(taggerhit);

        // prefilter events
        // check if MM within defined window
        if (!ParticleTypeDatabase::Proton.GetWindow(settings.mm_window_size).Round().Contains(comb.MissingMass))
            continue;
        h.steps->Fill("MM prefilter", 1);
        // try to reject unphysical proton candidates (proton < ~100MeV doesn't leave the target)
        // but only in BaF2 since the PbWO4 are a bit unreliable...
        if (taps->IsBaF2(comb.Proton->Candidate->FindCaloCluster()->CentralElement) && comb.Proton->Ek() < 50.)
            continue;
        h.steps->Fill("Sane p_E BaF2", 1);
        // reject photon candidates which have a too low energy
        // (photons > ~100MeV, lower likely split-off / noise / background)
        for (auto p : comb.Photons)
            if (p->Candidate->VetoEnergy < .3 && p->Ek() < 60.)
                continue;
        h.steps->Fill("#gamma E > 60", 1);
        // tighter PID timing
        for (auto p : comb.Photons)
            if (p->Candidate->VetoEnergy > .3 &&
                    (p->Candidate->FindVetoCluster()->Time < -10. || p->Candidate->FindVetoCluster()->Time > 32))
                continue;
        h.steps->Fill("PID timing", 1);
        // check if there are at least 2 PID entries
        if (std::count_if(comb.Photons.begin(), comb.Photons.end(),
                          [](TParticlePtr g){ return g->Candidate->VetoEnergy > .3; }) < 2)
            continue;
        h.steps->Fill("2 PIDs prefilter", 1);
        //test end

        // find best combination for each Tagger hit
        best_prob_fit = -std_ext::inf;
        // #combinations: binomial coefficient (n\\k)
        vector<double> IM_2g(3, std_ext::NaN);
        vector<double> IM_2g_fit(3, std_ext::NaN);

        //for (const auto& cand : selection) {
        for (const auto& cand : {comb}) {
            // do the fitting and check if the combination is better than the previous best
            if (!doFit_checkProb(taggerhit, cand, h, sig, best_prob_fit))
                continue;

            sig.DiscardedEk = cand.DiscardedEk;

            // run a kinematic fit with lepton candidates treated as charged pions and set the probability in the to-be-written tree
            sig.prob_antiPionFit = anti_pion_fit(taggerhit, cand);
            h.antiPionProb->Fill(sig.prob_antiPionFit);

            if (settings.less_plots())
                continue;
            utils::ParticleTools::FillIMCombinations(IM_2g.begin(), 2, cand.Photons);
            utils::ParticleTools::FillIMCombinations(IM_2g_fit.begin(), 2, sig.photons_kinfitted());
        }

        // only fill tree if a valid combination for the current Tagger hit was found
        if (!isfinite(best_prob_fit))
            continue;

        if (!settings.less_plots()) {
            for (const auto& im : IM_2g)
                h_subIM_2g->Fill(im, sig.TaggW);
            for (const auto& im : IM_2g_fit)
                h_subIM_2g_fit->Fill(im, sig.TaggW);
        }

        sig.Tree->Fill();
        h.steps->Fill("Tree filled", 1);
    }

    if (!isfinite(best_prob_fit))
        return;
    h.steps->Fill("best comb", 1);

    h_counts->Fill(chan_id.decaystring.c_str(), 1);

    if (settings.less_plots())
        return;

    auto get_veto_energies = [] (vector<TSimpleParticle> particles)
    {
        vector<double> veto_energies;
        for (const auto& p : particles)
            veto_energies.emplace_back(p.VetoE);

        return veto_energies;
    };

    const auto etap_fs = sig.photons();
    const auto veto_energies = get_veto_energies(etap_fs);

    // get sorted indices of the eta' final state according to their Veto energies
    const auto sorted_idx = std_ext::get_sorted_indices_desc(veto_energies);

    // do an anti pi0 cut on the combinations e+g and e-g
    // (assuming the photon deposited the least energy in the PIDs)
    if (Cuts_t::ANTI_PI0_CUT) {
        const interval<double> pion_cut(Cuts_t::ANTI_PI0_LOW, Cuts_t::ANTI_PI0_HIGH);
        LorentzVec pi0;
        const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};
        for (const auto pi0_comb : pi0_combs) {
            pi0 = LorentzVec({0,0,0,0});
            for (const auto idx : pi0_comb)
                pi0 += TParticle(ParticleTypeDatabase::Photon, etap_fs.at(sorted_idx[idx]));
            // apply an anti pion cut
            if (pion_cut.Contains(pi0.M()))
                return;
        }
        h.steps->Fill("anti #pi^{0} cut", 1);
    }

    const TSimpleParticle proton = sig.p;
    LorentzVec etap({0,0,0,0});
    for (const auto& g : etap_fs)
        etap += TParticle(ParticleTypeDatabase::Photon, g);
    h_protonVeto->Fill(proton.VetoE);
    h_pTheta->Fill(std_ext::radian_to_degree(proton.Theta()));

    ///\todo: still needed here? check final selection via ProtonPhotonCombs
    const auto idx1 = sorted_idx[0];
    const auto idx2 = sorted_idx[1];
    const auto l1 = etap_fs.at(idx1);
    const auto l2 = etap_fs.at(idx2);
    // suppress conversion decays
    if (sig.photons_vetoChannel().at(idx1) == sig.photons_vetoChannel().at(idx2))
        return;
    h.steps->Fill("distinct PID", 1);
    const double eeIM = (TParticle(ParticleTypeDatabase::eMinus, l1)
                         + TParticle(ParticleTypeDatabase::eMinus, l2)).M();
    h_IM2d->Fill(etap.M(), eeIM);
    h.IM2d->Fill(etap.M(), eeIM);

    // test effective cluster radius to distinguish between leptons and charged pions
    double effective_radius = sig.photons_effect_radius().at(idx1);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l1.Energy(), effective_radius);
    }
    effective_radius = sig.photons_effect_radius().at(idx2);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l2.Energy(), effective_radius);
    }
    double lateral_moment = sig.photons_lat_moment().at(idx1);
    if (isfinite(lateral_moment)) {
        h.lat_moment->Fill(lateral_moment);
        h.lat_moment_E->Fill(l1.Energy(), lateral_moment);
    }
    lateral_moment = sig.photons_lat_moment().at(idx2);
    if (isfinite(lateral_moment)) {
        h.lat_moment->Fill(lateral_moment);
        h.lat_moment_E->Fill(l2.Energy(), lateral_moment);
    }

    // test cluster size compared to energy
    h.cluster_size->Fill(l1.ClusterSize);
    h.cluster_size->Fill(l2.ClusterSize);
    h.cluster_size_E->Fill(l1.Energy(), l1.ClusterSize);
    h.cluster_size_E->Fill(l2.Energy(), l2.ClusterSize);

    h.etapIM_final->Fill(etap.M());
    h_etapIM_final->Fill(etap.M());
}

void EtapDalitz::ShowResult()
{
    for (auto& entry : channels)
        entry.second.Show();

    if (settings.less_plots())
        return;

    canvas(GetName()) << drawoption("colz") << h_IM2d << endc;

//    list<TH1*> hists;
//    for (auto& entry : channels) {
//        hists.push_back(entry.second.proton_E_theta);
//    }

//    hists.sort([](const TH1* a, const TH1* b) {return a->GetEntries() > b->GetEntries();});

//    int i=0;
//    for (auto& h : hists) {
//        c << h;
//        i++;
//        if (i>=9)
//            break;
//    }

//    c << endc;
}

bool EtapDalitz::doFit_checkProb(const TTaggerHit& taggerhit,
                                 const particle_comb_t& comb,
                                 PerChannel_t& h,
                                 SigTree_t& t,
                                 double& best_prob_fit)
{
    LorentzVec etap_kinfit({0,0,0,0});
    LorentzVec etap_treefit({0,0,0,0});
    LorentzVec etap_kinfit_freeZ({0,0,0,0});
    LorentzVec etap_treefit_freeZ({0,0,0,0});

    LorentzVec etap = comb.PhotonSum;
    h.etapIM->Fill(etap.M(), t.TaggW);
    h.MM->Fill(comb.MissingMass, t.TaggW);

    /* check what proton energy is predicted based on 4 momentum conservation,
     * set energy to measured if smaller than 350MeV (no punch-through) */
    {
        const auto beam_target = taggerhit.GetPhotonBeam() + LorentzVec::AtRest(ParticleTypeDatabase::Proton.Mass());
        sig.p_predictedEnergy = (beam_target - comb.PhotonSum).E - ParticleTypeDatabase::Proton.Mass();

        // set fitter uncertainty models depending on predicted energy and thus expected punch-through
        utils::UncertaintyModelPtr model;
        if (sig.p_predictedEnergy < 350.)
            model = sig.MCtrue ? model_MC_protonMeasured : model_data_protonMeasured;
        else
            model = sig.MCtrue ? model_MC : model_data;
        kinfit.SetUncertaintyModel(model);
        kinfit_freeZ.SetUncertaintyModel(model);
        treefitter_etap.SetUncertaintyModel(model);
        treefitter_etap_freeZ.SetUncertaintyModel(model);
    }


    /* start with the kinematic fitting */

    // treefit
    APLCON::Result_t treefit_result;

    treefitter_etap.PrepareFits(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

    // works this way because only one combination needs to be fitted
    treefitter_etap.NextFit(treefit_result);

    if (settings.use_treefit) {
        if (treefit_result.Status != APLCON::Result_Status_t::Success)
            return false;
        h.steps->Fill("treefit", 1);
    }

    // treefit free Z vertex
    APLCON::Result_t treefit_freeZ_result;

    treefitter_etap_freeZ.PrepareFits(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

    treefitter_etap_freeZ.NextFit(treefit_freeZ_result);


    // kinfit
    auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

    if (!settings.use_treefit) {
        if (kinfit_result.Status != APLCON::Result_Status_t::Success)
            return false;
        h.steps->Fill("kinfit", 1);
    }

    // kinfit free Z vertex
    auto kinfit_freeZ_result = kinfit_freeZ.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);


    const double kinfit_prob = kinfit_result.Probability;
    const double treefit_prob = treefit_result.Probability;

    h.treefitProb->Fill(treefit_prob);
    h.treefit_freeZ_prob->Fill(treefit_freeZ_result.Probability);
    h.kinfitProb->Fill(kinfit_prob);
    h.kinfit_freeZ_prob->Fill(kinfit_freeZ_result.Probability);

    // determine which probability should be used to find the best candidate combination
    const double prob = settings.use_treefit ? treefit_prob : kinfit_prob;

    if (Cuts_t::PROBABILITY_CUT) {
        if (prob < Cuts_t::PROBABILITY)
            return false;
        h.steps->Fill("probability", 1);
    }

    if (!settings.less_plots()) {
        h_etap->Fill(etap.E - etap.M(), std_ext::radian_to_degree(etap.Theta()), t.TaggW);
        h_proton->Fill(comb.Proton->E - comb.Proton->M(), std_ext::radian_to_degree(comb.Proton->Theta()), t.TaggW);
    }

    // check if a better probability has been found
    if (!std_ext::copy_if_greater(best_prob_fit, prob))
        return false;


    /* Gather relevant fitter information and update branches */

    t.set_proton_information(comb.Proton);
    t.set_photon_information(comb.Photons, true);  // store lateral moment and effective cluster radius
    t.set_additional_photon_information(comb.Photons);

    t.p_effect_radius = effective_radius(comb.Proton->Candidate);
    t.p_lat_moment    = lat_moment(comb.Proton->Candidate);

    t.etap = etap;
    t.mm   = comb.MissingMass;
    t.copl = std_ext::radian_to_degree(abs(etap.Phi() - comb.Proton->Phi())) - 180.;

    // now handle the different fitted particle information separately

    // kinfit
    if (kinfit_result.Status == APLCON::Result_Status_t::Success) {
        assert(kinfit.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : kinfit.GetFittedPhotons())
            etap_kinfit += *g;

        h.kinfit_ZVertex->Fill(kinfit.GetFittedZVertex());

        // update tree branches
        t.set_kinfit_information(kinfit, kinfit_result);
        t.etap_kinfit = etap_kinfit;
    }

    // treefit
    if (treefit_result.Status == APLCON::Result_Status_t::Success) {
        assert(treefitter_etap.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : treefitter_etap.GetFittedPhotons())
            etap_treefit += *g;

        h.treefit_ZVertex->Fill(treefitter_etap.GetFittedZVertex());

        // update tree branches
        t.set_treefit_information(treefitter_etap, treefit_result);
        t.etap_treefit = etap_treefit;
    }

    // handling free Z vertex fits starting here
    if (kinfit_freeZ_result.Status != APLCON::Result_Status_t::Success
            && treefit_freeZ_result.Status != APLCON::Result_Status_t::Success)
        return true;

    // update tree branches
    t.set_fit_freeZ_results(kinfit_freeZ, treefitter_etap_freeZ,
                            kinfit_freeZ_result, treefit_freeZ_result);

    // kinfit with free Z vertex
    if (kinfit_freeZ_result.Status == APLCON::Result_Status_t::Success) {
        assert(kinfit_freeZ.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : kinfit_freeZ.GetFittedPhotons())
            etap_kinfit_freeZ += *g;

        h.kinfit_freeZ_ZVertex->Fill(kinfit_freeZ.GetFittedZVertex());

        t.etap_kinfit_freeZ = etap_kinfit_freeZ;
    }

    // treefit with free Z vertex
    if (treefit_freeZ_result.Status == APLCON::Result_Status_t::Success) {
        assert(treefitter_etap_freeZ.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : treefitter_etap_freeZ.GetFittedPhotons())
            etap_treefit_freeZ += *g;

        h.treefit_freeZ_ZVertex->Fill(treefitter_etap_freeZ.GetFittedZVertex());

        t.etap_treefit_freeZ = etap_treefit_freeZ;
    }

    return true;
}

double EtapDalitz::anti_pion_fit(const TTaggerHit& taggerhit, const particle_comb_t& comb)
{
    fake_comb_t cand;
    cand.reset();
    cand.Proton = comb.Proton;

    auto leptons = get_sorted_indices_vetoE(comb.Photons);

    assert(leptons.size() == comb.Photons.size());

    // test the hypothesis of the two possible photon clusters with the highest veto energy to be charged pions
    cand.Photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::PiCharged, comb.Photons.at(leptons[0])->Candidate));
    cand.Photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::PiCharged, comb.Photons.at(leptons[1])->Candidate));
    cand.Photons.emplace_back(comb.Photons.at(leptons[2]));

    cand.calc_values(taggerhit);

    auto anti_fit_result = kinfit.DoFit(taggerhit.PhotonEnergy, cand.Proton, cand.Photons);

    if (anti_fit_result.Status != APLCON::Result_Status_t::Success)
        return -std_ext::inf;

    return anti_fit_result.Probability;
}

void EtapDalitzTools::count_clusters(const TCandidateList& cands, size_t& nCB, size_t& nTAPS)
{
    for (auto p : cands.get_iter())
        if (p->Detector & Detector_t::Type_t::CB)
            nCB++;
        else if (p->Detector & Detector_t::Type_t::TAPS)
            nTAPS++;
}

void EtapDalitzTools::channel_id(const TEvent& event, channel_id_t& chan_id)
{
    // assume data by default
    chan_id.production = "data";
    chan_id.decaystring = "data";
    chan_id.decay_name = "data";

    // get MC true channel information
    if (event.Reconstructed().ID.isSet(TID::Flags_t::MC)) {
        chan_id.production = std_ext::string_sanitize(utils::ParticleTools::GetProductionChannelString(event.MCTrue().ParticleTree).c_str());
        std_ext::remove_chars(chan_id.production, {'#', '{', '}', '^'});
        chan_id.decaystring = std_ext::string_sanitize(utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree).c_str());
        chan_id.decay_name = chan_id.decaystring;
        std_ext::remove_chars(chan_id.decay_name, {'#', '{', '}', '^'});
    }
}

EtapDalitz::PerChannel_t EtapDalitz::manage_channel_histograms_get_current(const bool MC, const TEvent& event)
{
    // check if the current production is known already, create new HistogramFactory otherwise
    auto prod = productions.find(chan_id.production);
    if (prod == productions.end()) {
        auto hf = new HistogramFactory(chan_id.production, HistFac, "");
        productions.insert({chan_id.production, *hf});
    }
    prod = productions.find(chan_id.production);
    auto hf = prod->second;

    // check if the decay channel is known already, if not insert it
    auto c = channels.find(chan_id.decaystring);
    if (c == channels.end())
        channels.insert({chan_id.decaystring, PerChannel_t(chan_id.decay_name, chan_id.decaystring, hf)});

    c = channels.find(chan_id.decaystring);
    if (MC && !Settings_t::get().less_plots())
        c->second.Fill(event.MCTrue());

    // return the histogram struct for the current channel
    return c->second;
}

bool EtapDalitzTools::q2_preselection(const TEventData& data, const double threshold = 50.) const
{
    auto mctrue_particles = utils::ParticleTypeList::Make(data.ParticleTree);
    auto particles = mctrue_particles.GetAll();
    // apply preselection condition only on channels with two leptons
    if (std::count_if(particles.begin(), particles.end(), [](TParticlePtr p){
                      return p->Type() == ParticleTypeDatabase::eCharged; }) != 2)
        return true;

    LorentzVec q2;
    for (auto p : mctrue_particles.GetAll())
        if (p->Type() == ParticleTypeDatabase::eCharged)
            q2 += *p;
    if (q2.M() > threshold)
        return true;

    return false;
}

void EtapDalitz::proton_tree::set_proton_information(const TParticlePtr proton)
{
    p             = TSimpleParticle(*proton);
    p_PSAangle    = proton->Candidate->FindCaloCluster()->GetPSAAngle();
    p_PSAradius   = proton->Candidate->FindCaloCluster()->GetPSARadius();
    p_detector    = getDetectorAsInt(proton->Candidate->Detector);
    p_centralElem = proton->Candidate->FindCaloCluster()->CentralElement;
    p_vetoChannel = 0;
    if (proton->Candidate->VetoEnergy) {
        p_vetoChannel = proton->Candidate->FindVetoCluster()->CentralElement;
        p_vetoTime    = proton->Candidate->FindVetoCluster()->Time;
    }
}

template <size_t N>
void EtapDalitz::photon_tree<N>::set_photon_information(const TParticleList& photons, const bool shower_shape)
{
    assert(photons.size() == N);
    this->photons() = TSimpleParticle::TransformParticleList(photons);

    for (size_t i = 0; i < photons.size(); ++i) {
        if (shower_shape) {
            photons_effect_radius().at(i) = effective_radius(photons.at(i)->Candidate);
            photons_lat_moment().at(i)    = lat_moment(photons.at(i)->Candidate);
        }
        photons_detector().at(i)      = getDetectorAsInt(photons.at(i)->Candidate->Detector);
        photons_centralElem().at(i)   = photons.at(i)->Candidate->FindCaloCluster()->CentralElement;
        photons_vetoChannel().at(i)   = 0;
        if (photons.at(i)->Candidate->VetoEnergy) {
            photons_vetoChannel().at(i) = photons.at(i)->Candidate->FindVetoCluster()->CentralElement;
            photons_vetoTime().at(i)    = photons.at(i)->Candidate->FindVetoCluster()->Time;
        }
    }
}

template<size_t Nphotons>
void EtapDalitz::fit_tree<Nphotons>::set_kinfit_information(const analysis::utils::KinFitter& kinfit,
                                                            const APLCON::Result_t& result)
{
    const auto kinfitted_photons = kinfit.GetFittedPhotons();
    const auto kinfit_particles  = kinfit.GetFitParticles();

    assert(kinfitted_photons.size() == Nphotons);

    kinfit_chi2        = result.ChiSquare;
    kinfit_probability = result.Probability;
    kinfit_iterations  = result.NIterations;
    kinfit_DoF         = result.NDoF;

    beam_E_kinfitted    = kinfit.GetFittedBeamE();
    beam_kinfit_E_pull  = kinfit.GetBeamEPull();
    kinfit_ZVertex      = kinfit.GetFittedZVertex();
    kinfit_ZVertex_pull = kinfit.GetZVertexPull();

    p_kinfitted = *kinfit.GetFittedProton();

    p_kinfit_E_pull     = kinfit_particles.at(0).GetPulls().at(0);
    p_kinfit_theta_pull = kinfit_particles.at(0).GetPulls().at(1);
    p_kinfit_phi_pull   = kinfit_particles.at(0).GetPulls().at(2);

    for (size_t i = 0; i < Nphotons; ++i) {
        photons_kinfitted().at(i) = *(kinfitted_photons.at(i));

        photon_kinfit_E_pulls().at(i)     = kinfit_particles.at(i+1).GetPulls().at(0);
        photon_kinfit_theta_pulls().at(i) = kinfit_particles.at(i+1).GetPulls().at(1);
        photon_kinfit_phi_pulls().at(i)   = kinfit_particles.at(i+1).GetPulls().at(2);
    }
}

template<size_t Nphotons>
void EtapDalitz::fit_tree<Nphotons>::set_treefit_information(const analysis::utils::TreeFitter& treefit,
                                                             const APLCON::Result_t& result)
{
    const auto treefitted_photons = treefit.GetFittedPhotons();
    const auto treefit_particles  = treefit.GetFitParticles();

    assert(treefitted_photons.size() == Nphotons);

    treefit_chi2        = result.ChiSquare;
    treefit_probability = result.Probability;
    treefit_iterations  = result.NIterations;
    treefit_DoF         = result.NDoF;

    beam_E_treefitted    = treefit.GetFittedBeamE();
    beam_treefit_E_pull  = treefit.GetBeamEPull();
    treefit_ZVertex      = treefit.GetFittedZVertex();
    treefit_ZVertex_pull = treefit.GetZVertexPull();

    p_treefitted = *treefit.GetFittedProton();

    p_treefit_E_pull     = treefit_particles.at(0).GetPulls().at(0);
    p_treefit_theta_pull = treefit_particles.at(0).GetPulls().at(1);
    p_treefit_phi_pull   = treefit_particles.at(0).GetPulls().at(2);

    for (size_t i = 0; i < Nphotons; ++i) {
        photons_treefitted().at(i) = *(treefitted_photons.at(i));

        photon_treefit_E_pulls().at(i)     = treefit_particles.at(i+1).GetPulls().at(0);
        photon_treefit_theta_pulls().at(i) = treefit_particles.at(i+1).GetPulls().at(1);
        photon_treefit_phi_pulls().at(i)   = treefit_particles.at(i+1).GetPulls().at(2);
    }
}

template<size_t Nphotons>
void EtapDalitz::fit_freeZ_tree<Nphotons>::set_fit_freeZ_results(const analysis::utils::KinFitter& kinfit,
                                                                 const analysis::utils::TreeFitter& treefit,
                                                                 const APLCON::Result_t& kinfit_result,
                                                                 const APLCON::Result_t& treefit_result)
{
    // kinfit with free Z vertex
    if (kinfit_result.Status == APLCON::Result_Status_t::Success) {
        const auto kinfit_freeZ_photons   = kinfit.GetFittedPhotons();
        const auto kinfit_freeZ_particles = kinfit.GetFitParticles();

        assert(kinfit_freeZ_photons.size() == Nphotons);

        kinfit_freeZ_chi2        = kinfit_result.ChiSquare;
        kinfit_freeZ_probability = kinfit_result.Probability;
        kinfit_freeZ_iterations  = kinfit_result.NIterations;
        kinfit_freeZ_DoF         = kinfit_result.NDoF;

        beam_E_kinfit_freeZ       = kinfit.GetFittedBeamE();
        beam_kinfit_freeZ_E_pull  = kinfit.GetBeamEPull();
        kinfit_freeZ_ZVertex      = kinfit.GetFittedZVertex();
        kinfit_freeZ_ZVertex_pull = kinfit.GetZVertexPull();

        p_kinfit_freeZ = *kinfit.GetFittedProton();

        p_kinfit_freeZ_theta_pull = kinfit_freeZ_particles.at(0).GetPulls().at(1);
        p_kinfit_freeZ_phi_pull   = kinfit_freeZ_particles.at(0).GetPulls().at(2);

        for (size_t i = 0; i < Nphotons; ++i) {
            photons_kinfit_freeZ().at(i) = *(kinfit_freeZ_photons.at(i));

            photon_kinfit_freeZ_E_pulls().at(i)     = kinfit_freeZ_particles.at(i+1).GetPulls().at(0);
            photon_kinfit_freeZ_theta_pulls().at(i) = kinfit_freeZ_particles.at(i+1).GetPulls().at(1);
            photon_kinfit_freeZ_phi_pulls().at(i)   = kinfit_freeZ_particles.at(i+1).GetPulls().at(2);
        }
    }

    // treefit with free Z vertex
    if (treefit_result.Status == APLCON::Result_Status_t::Success) {
        const auto treefit_freeZ_photons   = treefit.GetFittedPhotons();
        const auto treefit_freeZ_particles = treefit.GetFitParticles();

        assert(treefit_freeZ_photons.size() == Nphotons);

        treefit_freeZ_chi2        = treefit_result.ChiSquare;
        treefit_freeZ_probability = treefit_result.Probability;
        treefit_freeZ_iterations  = treefit_result.NIterations;
        treefit_freeZ_DoF         = treefit_result.NDoF;

        beam_E_treefit_freeZ       = treefit.GetFittedBeamE();
        beam_treefit_freeZ_E_pull  = treefit.GetBeamEPull();
        treefit_freeZ_ZVertex      = treefit.GetFittedZVertex();
        treefit_freeZ_ZVertex_pull = treefit.GetZVertexPull();

        p_treefit_freeZ = *treefit.GetFittedProton();

        p_treefit_freeZ_theta_pull = treefit_freeZ_particles.at(0).GetPulls().at(1);
        p_treefit_freeZ_phi_pull   = treefit_freeZ_particles.at(0).GetPulls().at(2);

        for (size_t i = 0; i < Nphotons; ++i) {
            photons_treefit_freeZ().at(i) = *(treefit_freeZ_photons.at(i));

            photon_treefit_freeZ_E_pulls().at(i)     = treefit_freeZ_particles.at(i+1).GetPulls().at(0);
            photon_treefit_freeZ_theta_pulls().at(i) = treefit_freeZ_particles.at(i+1).GetPulls().at(1);
            photon_treefit_freeZ_phi_pulls().at(i)   = treefit_freeZ_particles.at(i+1).GetPulls().at(2);
        }
    }
}

void EtapDalitz::SigTree_t::set_additional_photon_information(const TParticleList& photons)
{
    this->photons() = TSimpleParticle::TransformParticleList(photons);
    for (size_t i = 0; i < photons.size(); ++i) {
        photons_PSAangle().at(i)  = photons.at(i)->Candidate->FindCaloCluster()->GetPSAAngle();
        photons_PSAradius().at(i) = photons.at(i)->Candidate->FindCaloCluster()->GetPSARadius();
    }
}

void EtapDalitz::common_tree::init()
{
    // init values which are constant the whole event
    nCands    = 0;
    channel   = 0;

    CBSumE    = -std_ext::inf;
    CBAvgTime = std_ext::NaN;

    TaggW     = std_ext::NaN;
    TaggE     = std_ext::NaN;
    TaggT     = std_ext::NaN;
    TaggCh    = 0;
}

void EtapDalitz::common_tree::reset()
{
    // only reset values which change during an event
    TaggW  = std_ext::NaN;
    TaggE  = std_ext::NaN;
    TaggT  = std_ext::NaN;
    TaggCh = 0;
}

void EtapDalitz::proton_tree::reset()
{
    p             = TSimpleParticle();
    p_PSAangle    = std_ext::NaN;
    p_PSAradius   = std_ext::NaN;
    p_detector    = -1;
    p_centralElem = 0;
    p_vetoChannel = 0;
    p_vetoTime    = -std_ext::inf;
}

template <size_t N>
void EtapDalitz::photon_tree<N>::reset()
{
    for (size_t i = 0; i < N; ++i) {
        photons_effect_radius().at(i) = -std_ext::inf;
        photons_lat_moment().at(i)    = -std_ext::inf;
        photons_detector().at(i)      = -1;
        photons_centralElem().at(i)   = 0;
        photons_vetoChannel().at(i)   = 0;
        photons_vetoTime().at(i)      = -std_ext::inf;
    }
}

template <size_t Nphotons>
void EtapDalitz::fit_tree<Nphotons>::reset()
{
    beam_E_kinfitted    = -std_ext::inf;
    beam_E_treefitted   = -std_ext::inf;
    kinfit_chi2         = std_ext::NaN;
    kinfit_probability  = std_ext::NaN;
    treefit_chi2        = std_ext::NaN;
    treefit_probability = std_ext::NaN;

    etap_kinfit  = LorentzVec();
    etap_treefit = LorentzVec();
}

template <size_t Nphotons>
void EtapDalitz::fit_freeZ_tree<Nphotons>::reset()
{
    beam_E_kinfit_freeZ       = -std_ext::inf;
    beam_E_treefit_freeZ      = -std_ext::inf;
    kinfit_freeZ_chi2         = std_ext::NaN;
    kinfit_freeZ_probability  = std_ext::NaN;
    treefit_freeZ_chi2        = std_ext::NaN;
    treefit_freeZ_probability = std_ext::NaN;

    etap_kinfit_freeZ  = LorentzVec();
    etap_treefit_freeZ = LorentzVec();
}

void EtapDalitz::SigTree_t::init()
{
    common_tree::init();
}

void EtapDalitz::SigTree_t::reset()
{
    common_tree::reset();
    proton_tree::reset();
    photon_tree::reset();
    fit_tree::reset();
    fit_freeZ_tree::reset();

    for (size_t i = 0; i < this->photons().size(); ++i) {
        photons_PSAangle().at(i)  = std_ext::NaN;
        photons_PSAradius().at(i) = std_ext::NaN;
    }
    prob_antiPionFit = std_ext::NaN;
    p_predictedEnergy = std_ext::NaN;
}

void EtapDalitz::RefTree_t::init()
{
    common_tree::init();
}

void EtapDalitz::RefTree_t::reset()
{
    common_tree::reset();
    proton_tree::reset();
    photon_tree::reset();
    fit_tree::reset();
}

EtapDalitz::ReactionChannel_t::ReactionChannel_t(const string &n):
    name(n)
{}

EtapDalitz::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<EtapDalitz::decaytree_t> &t, const short c):
    name(utils::ParticleTools::GetDecayString(t)),
    tree(t),
    color(c)
{}

EtapDalitz::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<EtapDalitz::decaytree_t> &t, const string &n, const short c):
    name(n),
    tree(t),
    color(c)
{}

EtapDalitz::ReactionChannelList_t EtapDalitz::makeChannels()
{
    ReactionChannelList_t m;

    m.channels[0] = ReactionChannel_t("Data");
    m.channels[1] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_eeg), kRed};  // signal
    m.channels[2] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g), kRed+1};  // reference
    m.channels[3] = ReactionChannel_t("Sum MC");
    m.channels[4] = ReactionChannel_t("MC BackG");

    m.channels[10] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),           "#pi^{0} #pi^{0} #rightarrow 4#gamma", kGreen-4};
    m.channels[11] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),           "#pi^{0} #eta #rightarrow 4#gamma", kAzure+1};
    m.channels[12] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gRho_gPiPi), "#eta' #rightarrow #rho #gamma", kOrange+6};
    m.channels[13] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g),         "#eta' #rightarrow #gamma #gamma", kOrange};
    m.channels[14] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_2ggEpEm),      "#pi^{0} #pi^{0} #rightarrow 2#gamma e^{+} e^{-} #gamma", kSpring+10};
    m.channels[15] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Rho_PiPi),            "#rho #rightarrow #pi^{+} #pi^{-}", kBlue+1};
    m.channels[16] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),              "#pi^{0} #rightarrow #gamma #gamma", kTeal};
    m.channels[17] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0PiPi_2gPiPi),      "#pi^{0} #pi^{+} #pi^{-}", kViolet-4};
    m.channels[18] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_2g),              "#eta #rightarrow #gamma #gamma", kOrange+4};
    m.channels[m.other_index] = ReactionChannel_t(nullptr, "Others", kCyan-6);

    return m;
}

unsigned EtapDalitz::ReactionChannelList_t::identify(const ant::TParticleTree_t& tree) const
{
    if (!tree)
        return 0;

    for (const auto& c : channels) {

        if (!c.second.tree)
            continue;

        if (tree->IsEqual(c.second.tree, utils::ParticleTools::MatchByParticleName))
            return c.first;
    }

    return other_index;
}

const EtapDalitz::ReactionChannelList_t EtapDalitz::reaction_channels = EtapDalitz::makeChannels();

const unsigned EtapDalitz::ReactionChannelList_t::other_index = 100;


/* Reference channel analysis */
Etap2g::Etap2g(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    model_data(utils::UncertaintyModels::Interpolated::makeAndLoad(
                   utils::UncertaintyModels::Interpolated::Type_t::Data,
                   // use Sergey as starting point
                   make_shared<utils::UncertaintyModels::FitterSergey>()
                   )),
    model_MC(utils::UncertaintyModels::Interpolated::makeAndLoad(
                 utils::UncertaintyModels::Interpolated::Type_t::MC,
                 // use Sergey as starting point
                 make_shared<utils::UncertaintyModels::FitterSergey>()
                 )),
    kinfit(nullptr,
           opts->HasOption("SigmaZ"), EtapDalitz::MakeFitSettings(20)
           ),
    treefitter_etap(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g),
                    nullptr, opts->HasOption("SigmaZ"), {}, EtapDalitz::MakeFitSettings(20)
                    )
{
    if (opts->HasOption("SigmaZ")) {
        double sigma_z = 0.;
        std_ext::copy_if_greater(sigma_z, opts->Get<double>("SigmaZ", 0.));
        kinfit.SetZVertexSigma(sigma_z);
        treefitter_etap.SetZVertexSigma(sigma_z);
    }
}

void Etap2g::ProcessEvent(const TEvent& event, manager_t&)
{
    Process(event);
}

void Etap2g::Process(const TEvent& event)
{
    triggersimu.ProcessEvent(event);

    const bool MC = event.Reconstructed().ID.isSet(TID::Flags_t::MC);
    const auto& cands = event.Reconstructed().Candidates;

    if (t->nCands != Cuts_t::N_FINAL_STATE)
        return;

    // set fitter uncertainty models
    kinfit.SetUncertaintyModel(MC ? model_MC : model_data);
    treefitter_etap.SetUncertaintyModel(MC ? model_MC : model_data);

    /* test combinations to determine best proton candidate */
    if (TEST_COMBS) {
        utils::ProtonPhotonCombs proton_photons(cands);
        double best_prob_fit = -std_ext::inf;

        // loop over all tagger hits
        for (const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {
            promptrandom->SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
            if (promptrandom->State() == PromptRandom::Case::Outside)
                continue;

            t->TaggW = promptrandom->FillWeight();
            t->TaggE = taggerhit.PhotonEnergy;
            t->TaggT = taggerhit.Time;
            t->TaggCh = taggerhit.Channel;

            // find best combination for each Tagger hit
            best_prob_fit = -std_ext::inf;

            particle_combs_t selection = proton_photons()
                    .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(Cuts_t::MM_WINDOW_SIZE).Round())  // MM cut on proton mass
                    .FilterCustom([=] (const particle_comb_t& p) {
                        // ensure the possible proton candidate is kinematically allowed
                        if (std_ext::radian_to_degree(p.Proton->Theta()) > Cuts_t::MAX_PROTON_THETA)
                            return true;
                        return false;
                    }, "proton #vartheta");
//                    .FilterCustom([] (const particle_comb_t& p) {
//                        // require that PID entries do not have more than .3 MeV
//                        if (std::count_if(p.Photons.begin(), p.Photons.end(),
//                                          [](TParticlePtr g){ return g->Candidate->VetoEnergy < .3; }) < 2)
//                            return true;
//                        return false;
//                    }, "PID energy");

            if (selection.empty())
                continue;

            for (const auto& cand : selection) {
                // do the fitting and check if the combination is better than the previous best
                if (!doFit_checkProb(taggerhit, cand.Proton, cand.Photons, best_prob_fit))
                    continue;
            }

            // only fill tree if a valid combination for the current Tagger hit was found
            if (!isfinite(best_prob_fit))
                continue;

            t->Tree->Fill();
        }

        return;
    }


    /* use simple selection requiring proton in TAPS and 2 photons in CB */
    TParticlePtr proton;
    TParticleList photons;

    if (!simple2CB1TAPS(cands, proton, photons))
        return;

    for (const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {  // loop over all tagger hits
        t->reset();
        promptrandom->SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if (promptrandom->State() == PromptRandom::Case::Outside)
            continue;

        t->TaggW = promptrandom->FillWeight();
        t->TaggE = taggerhit.PhotonEnergy;
        t->TaggT = taggerhit.Time;
        t->TaggCh = taggerhit.Channel;

        /* kinematic fitting */
        // treefit
        APLCON::Result_t treefit_result;

        treefitter_etap.PrepareFits(taggerhit.PhotonEnergy, proton, photons);

        treefitter_etap.NextFit(treefit_result);  // no loop, just one possible combination

        if (USE_TREEFIT && treefit_result.Status != APLCON::Result_Status_t::Success)
            return;

        // kinfit

        auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, proton, photons);

        if (USE_KINFIT && kinfit_result.Status != APLCON::Result_Status_t::Success)
            return;

        // fill the tree with the fitted values
        fill_tree(treefit_result, kinfit_result, proton, photons);
        t->Tree->Fill();
    }

}

void Etap2g::fill_tree(const APLCON::Result_t& treefit_result,
                       const APLCON::Result_t& kinfit_result,
                       const TParticlePtr proton,
                       const TParticleList& photons)
{
    LorentzVec etap;
    LorentzVec etap_kinfit;
    LorentzVec etap_treefit;
    /* check if the fits converged and fill the trees accordingly */

    // update branches with general particle and fitter information
    etap = sumlv(photons.begin(), photons.end());

    t->set_proton_information(proton);
    t->set_photon_information(photons, false);  // don't store lateral moment and effective cluster radius
    t->etap = etap;

    // kinfit
    if (kinfit_result.Status == APLCON::Result_Status_t::Success) {
        assert(kinfit.GetFitParticles().size() == Cuts_t::N_FINAL_STATE);

        auto kinfit_photons = kinfit.GetFittedPhotons();

        etap_kinfit = sumlv(kinfit_photons.begin(), kinfit_photons.end());

        // update tree branches
        t->set_kinfit_information(kinfit, kinfit_result);
        t->etap_kinfit = etap_kinfit;
    }

    // treefit
    if (treefit_result.Status == APLCON::Result_Status_t::Success) {
        assert(treefitter_etap.GetFitParticles().size() == Cuts_t::N_FINAL_STATE);

        auto treefit_photons = treefitter_etap.GetFittedPhotons();

        etap_treefit = sumlv(treefit_photons.begin(), treefit_photons.end());

        // update tree branches
        t->set_treefit_information(treefitter_etap, treefit_result);
        t->etap_treefit = etap_treefit;
    }
}

bool Etap2g::simple2CB1TAPS(const TCandidateList& cands,
                            TParticlePtr& proton,
                            TParticleList& photons)
{
    size_t nCB = 0, nTAPS = 0;
    photons.clear();

    for (auto p : cands.get_iter())
        if (p->Detector & Detector_t::Type_t::CB && nCB++ < 2)
            photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, p));
        else if (p->Detector & Detector_t::Type_t::TAPS && nTAPS++ < 1)
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, p);
        else if (nCB > 2 || nTAPS > 1)
            return false;

    return true;
}

bool Etap2g::doFit_checkProb(const TTaggerHit& taggerhit,
                             const TParticlePtr proton,
                             const TParticleList& photons,
                             double& best_prob_fit)
{
    LorentzVec etap({0,0,0,0});

    for (const auto& g : photons)
        etap += *g;

    /* kinematical checks to reduce computing time */
    const interval<double> mm = ParticleTypeDatabase::Proton.GetWindow(Cuts_t::MM_WINDOW_SIZE);

    LorentzVec missing = taggerhit.GetPhotonBeam() + LorentzVec::AtRest(ParticleTypeDatabase::Proton.Mass());
    missing -= etap;
    if (!mm.Contains(missing.M()))
        return false;


    /* now start with the kinematic fitting */
    // treefit
    APLCON::Result_t treefit_result;

    treefitter_etap.PrepareFits(taggerhit.PhotonEnergy, proton, photons);

    // works this way because only one combination needs to be fitted
    while (treefitter_etap.NextFit(treefit_result))
        if (treefit_result.Status != APLCON::Result_Status_t::Success)
            continue;

    if (USE_TREEFIT)
        if (treefit_result.Status != APLCON::Result_Status_t::Success)
            return false;

    // kinfit

    auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, proton, photons);

    if (!USE_TREEFIT)
        if (kinfit_result.Status != APLCON::Result_Status_t::Success)
            return false;


    // determine which probability should be used to find the best candidate combination
    const double prob = USE_TREEFIT ? treefit_result.Probability : kinfit_result.Probability;

    if (Cuts_t::PROBABILITY_CUT) {
        if (prob < Cuts_t::PROBABILITY)
            return false;
    }

    if (!std_ext::copy_if_greater(best_prob_fit, prob))
        return false;

    fill_tree(treefit_result, kinfit_result, proton, photons);

    return true;
}

void Etap2g::setPromptRandom(PromptRandom::Switch& prs)
{
    promptrandom = &prs;
}

void Etap2g::linkTree(RefTree_t& ref)
{
    t = &ref;
}


AUTO_REGISTER_PHYSICS(EtapDalitz)
