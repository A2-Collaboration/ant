#include "etaprime_dalitz.h"

#include "utils/Combinatorics.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"


using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

template<typename T>
void shift_right(std::vector<T>& v)
{
    std::rotate(v.begin(), v.end() -1, v.end());
}

template<typename iter>
LorentzVec EtapDalitz::sumlv(iter start, iter end) {
    LorentzVec s;
    while (start != end) {
        s += **(start);
        ++start;
    }
    return s;
}

void EtapDalitz::remove_char(std::string& str, char ch)
{
    str.erase(std::remove(str.begin(), str.end(), ch), str.end());
}

void EtapDalitz::remove_chars(std::string& str, std::initializer_list<char> chars)
{
    for (const auto ch : chars)
        remove_char(str, ch);
}

double EtapDalitz::calc_effective_radius(const TCandidatePtr cand) const
{
    TClusterHitList crystals = cand->FindCaloCluster()->Hits;
    if (crystals.size() < 3)
        return std_ext::NaN;
    double effR = 0, e = 0;
    vec3 central = cb->GetPosition(cand->FindCaloCluster()->CentralElement);
    for (TClusterHit crystal : crystals) {
        const double r = std_ext::radian_to_degree(central.Angle(cb->GetPosition(crystal.Channel)));
        effR += r*r*crystal.Energy;
        e += crystal.Energy;
    }
    effR /= e;
    return sqrt(effR);
}

double EtapDalitz::lat_moment(const TCandidatePtr cand) const
{
    return clustertools.LateralMoment(*(cand->FindCaloCluster()));
}

ParticleTypeTree EtapDalitz::base_tree()
{
    ParticleTypeTree t = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(ParticleTypeDatabase::BeamProton);
    t->CreateDaughter(ParticleTypeDatabase::Proton);
    return t;
}

ParticleTypeTree EtapDalitz::etap_3g()
{
    auto t = base_tree();
    auto etap = t->CreateDaughter(ParticleTypeDatabase::EtaPrime);
    etap->CreateDaughter(ParticleTypeDatabase::Photon);
    etap->CreateDaughter(ParticleTypeDatabase::Photon);
    etap->CreateDaughter(ParticleTypeDatabase::Photon);
    return t;
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
    etapIM_kinfit = hf.makeTH1D(title + " IM #eta' kinfitted", "IM [MeV]", "#", energy, name + " etapIM_kinfit");
    etapIM_kinfit_freeZ = hf.makeTH1D(title + " IM #eta' free Z kinfitted", "IM [MeV]", "#", energy, name + " etapIM_kinfit_freeZ");
    etapIM_treefit = hf.makeTH1D(title + " IM #eta' treefitted", "IM [MeV]", "#", energy, name + " etapIM_treefit");
    etapIM_treefit_freeZ = hf.makeTH1D(title + " IM #eta' free Z treefitted", "IM [MeV]", "#", energy, name + " etapIM_treefit_freeZ");
    etapIM_cand = hf.makeTH1D(title + " IM #eta' candidates", "IM [MeV]", "#", energy, name + " etapIM_cand");
    etapIM_final = hf.makeTH1D(title + " IM #eta' final", "IM [MeV]", "#", energy, name + " etapIM_final");
    IM2d = hf.makeTH2D(title + " IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), name + " IM2d");
    MM = hf.makeTH1D(title + " Missing Mass proton", "MM [MeV]", "#", BinSettings(1600), name + " MM");
    hCopl = hf.makeTH1D(title + " Coplanarity #eta - proton all comb", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), name + " hCopl");
    hCopl_final = hf.makeTH1D(title + " Coplanarity #eta - proton final", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), name + " hCopl_final");
    treefitChi2 = hf.makeTH1D(title + " treefitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " treefitChi2");
    treefitProb = hf.makeTH1D(title + " treefitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " treefitProb");
    treefitIter = hf.makeTH1D(title + " treefitted # Iterations", "#iterations", "#", BinSettings(20), name + " treefitIter");
    treefit_freeZ_chi2 = hf.makeTH1D(title + " free Z treefitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " treefit_freeZ_chi2");
    treefit_freeZ_prob = hf.makeTH1D(title + " free Z treefitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " treefit_freeZ_prob");
    treefit_freeZ_iter = hf.makeTH1D(title + " free Z treefitted # Iterations", "#iterations", "#", BinSettings(20), name + " treefit_freeZ_iter");
    kinfitChi2 = hf.makeTH1D(title + " kinfitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " kinfitChi2");
    kinfitProb = hf.makeTH1D(title + " kinfitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " kinfitProb");
    kinfitIter = hf.makeTH1D(title + " kinfitted # Iterations", "#iterations", "#", BinSettings(20), name + " kinfitIter");
    kinfit_freeZ_chi2 = hf.makeTH1D(title + " free Z kinfitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " kinfit_freeZ_chi2");
    kinfit_freeZ_prob = hf.makeTH1D(title + " free Z kinfitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " kinfit_freeZ_prob");
    kinfit_freeZ_iter = hf.makeTH1D(title + " free Z kinfitted # Iterations", "#iterations", "#", BinSettings(20), name + " kinfit_freeZ_iter");
    kinfit_ZVertex = hf.makeTH1D(title + " kinfitted Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " kinfit_ZVertex");
    kinfit_freeZ_ZVertex = hf.makeTH1D(title + " free Z kinfitted Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " kinfit_freeZ_ZVertex");
    treefit_ZVertex = hf.makeTH1D(title + " treefitted Z Vertex", "z [cm]", "#", BinSettings(300, -30, 30), name + " treefit_ZVertex");
    treefit_freeZ_ZVertex = hf.makeTH1D(title + " free Z treefitted Z Vertex", "z [cm]", "#", BinSettings(300, -30, 30), name + " treefit_freeZ_ZVertex");
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
    //canvas("Per Channel: " + title) << drawoption("colz") << IM2d << endc;
    canvas("Per Channel: " + title) << steps
                                    << etapIM_kinfit
                                    << etapIM_treefit
                                    << etapIM_final
                                    << hCopl_final
                                    << kinfitChi2
                                    << kinfitProb
                                    << kinfit_freeZ_chi2
                                    << kinfit_freeZ_prob
                                    << treefitChi2
                                    << treefitProb
                                    << endc;
}

void EtapDalitz::PerChannel_t::Fill(const TEventData& d)
{
    auto particles = d.ParticleTree ?
                         utils::ParticleTypeList::Make(d.ParticleTree) :
                         utils::ParticleTypeList::Make(d.Candidates);
    const auto& protons = particles.Get(ParticleTypeDatabase::Proton);
    if (!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), p->Theta()*TMath::RadToDeg());
    }
}

static TVector2 getPSAVector(const TParticlePtr& p)
{
    if (p->Candidate) {
        const auto cluster = p->Candidate->FindCaloCluster();
        if (cluster)
            return {cluster->Energy, cluster->ShortEnergy};
    }

    throw std::runtime_error("Incomplete Particle without candiate or CaloCluster");
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
    reference(opts->Get<bool>("reference", 0)),
    reference_only(opts->Get<bool>("reference_only", 0)),
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

    cb = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    sig.CreateBranches(HistFac.makeTTree("signal"));
    if (reference || reference_only) {
        if (reference)
            LOG(INFO) << "Reference channel included in analysis";
        if (reference_only)
            LOG(INFO) << "Only Reference channel will be analysed";
        ref.CreateBranches(HistFac.makeTTree("ref"));
        etap2g = new Etap2g("Etap2g", opts);
        etap2g->linkTree(ref);
        etap2g->setPromptRandom(promptrandom);
    }

    const BinSettings tagger_time_bins(2000, -200, 200);

    h_tagger_time = HistFac.makeTH1D("Tagger Time", "t [ns]", "#", tagger_time_bins, "h_tagger_time");
    h_tagger_time_CBavg = HistFac.makeTH1D("Tagger Time - CB avg time", "t [ns]", "#", tagger_time_bins, "h_tagger_time_CBavg");

    h_counts = HistFac.makeTH1D("Events per Channel", "channel", "#", BinSettings(20), "h_counts");
    h_nCands = HistFac.makeTH1D("Number of Candidates", "#Candidates", "#", BinSettings(30), "h_nCands");
    h_cluster_CB = HistFac.makeTH1D("#Cluster CB", "#Cluster", "#", BinSettings(20), "h_cluster_CB");
    h_cluster_TAPS = HistFac.makeTH1D("#Cluster TAPS", "#Cluster", "#", BinSettings(20), "h_cluster_TAPS");
    missed_channels = HistFac.makeTH1D("Unlisted Channels", "", "Total Events seen", BinSettings(20), "missed_channels");
    found_channels  = HistFac.makeTH1D("Listed Channels", "", "Total Events seen", BinSettings(20), "found_channels");

    const BinSettings energybins(1000, 0, 10);

    h_pTheta = HistFac.makeTH1D("#vartheta proton candidate", "#vartheta_{p} [#circ]", "#", BinSettings(720, 0, 180), "h_pTheta");
    h_protonVeto = HistFac.makeTH1D("Veto energy identified proton", "Veto [MeV]", "#", energybins, "h_protonVeto");
    h_etapIM_final = HistFac.makeTH1D("IM #eta' final", "IM [MeV]", "#", BinSettings(1200), "h_etapIM_final");
    h_IM2d = HistFac.makeTH2D("IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), "h_IM2d");
    h_etap = HistFac.makeTH2D("Kinematics #eta'", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(360, 0, 180), "h_etap");
    h_proton = HistFac.makeTH2D("Kinematics p", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(160, 0, 80), "h_proton");

    // set sigma to 0 for unmeasured --> free z vertex
    kinfit_freeZ.SetZVertexSigma(0);
    treefitter_etap_freeZ.SetZVertexSigma(0);

    if (opts->HasOption("SigmaZ")) {
        double sigma_z = 0.;
        std_ext::copy_if_greater(sigma_z, opts->Get<double>("SigmaZ", 0.));
        LOG(INFO) << "Fit Z vertex enabled with sigma = " << sigma_z;
        kinfit.SetZVertexSigma(sigma_z);
        treefitter_etap.SetZVertexSigma(sigma_z);
    }
}

void EtapDalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

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

    std::string production = "data";
    std::string decaystring = "data";
    std::string decay_name = "data";
    if (MC) {
        production = std_ext::string_sanitize(utils::ParticleTools::GetProductionChannelString(event.MCTrue().ParticleTree).c_str());
        remove_chars(production, {'#', '{', '}', '^'});
        decaystring = std_ext::string_sanitize(utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree).c_str());
        decay_name = decaystring;
        remove_chars(decay_name, {'#', '{', '}', '^'});
    }

    auto prod = productions.find(production);
    if (prod == productions.end()) {
        auto hf = new HistogramFactory(production, HistFac, "");
        productions.insert({production, *hf});
    }
    prod = productions.find(production);
    auto hf = prod->second;

    auto c = channels.find(decaystring);
    if (c == channels.end())
        channels.insert({decaystring, PerChannel_t(decay_name, decaystring, hf)});

    c = channels.find(decaystring);
    if (MC)
        c->second.Fill(event.MCTrue());
    auto h = c->second;

    const auto& cands = data.Candidates;
    //const auto nCandidates = cands.size();
    sig.nCands = cands.size();
    h_nCands->Fill(sig.nCands);
    h.steps->Fill("seen", 1);

    // histogram amount of CB and TAPS clusters
    count_clusters(cands);

    if(!triggersimu.HasTriggered())
        return;

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

    if (reference || reference_only) {
        ref.MCtrue = sig.MCtrue;
        ref.channel = sig.channel;
        ref.trueZVertex = sig.trueZVertex;
        ref.nCands = sig.nCands;
        ref.CBSumE = sig.CBSumE;
        ref.CBAvgTime = sig.CBAvgTime;

        etap2g->Process(event);

        if (reference_only)
            return;
    }


    if (cands.size() != N_FINAL_STATE)
        return;
    h.steps->Fill("#cands", 1);

    // q2 preselection on MC data
    if (MC && Q2_PRESELECTION) {
        if (!q2_preselection(event.MCTrue(), 50.))
            return;
        h.steps->Fill("MC q2 > 50", 1);
    }

    TLorentzVector etap;
    TParticlePtr proton;
    //const interval<double> etap_im({ETAP_IM-ETAP_SIGMA, ETAP_IM+ETAP_SIGMA});
    TCandidatePtrList comb;
    for (auto p : cands.get_iter())
        comb.emplace_back(p);


    // require at least 2 candidates with PID/Veto entries
    if (std::count_if(comb.begin(), comb.end(), [](TCandidatePtr c){ return c->VetoEnergy; }) < 2)
        return;
    h.steps->Fill("#Veto", 1);

    TParticleList photons;
    double best_prob_fit = -std_ext::inf;
    size_t best_comb_fit = cands.size();
    for (const TTaggerHit& taggerhit : data.TaggerHits) {  // loop over all tagger hits
        if (!MC) {
            h_tagger_time->Fill(taggerhit.Time);
            h_tagger_time_CBavg->Fill(taggerhit.Time - sig.CBAvgTime);
        }

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h.steps->Fill("time window", 1);

        sig.TaggW = promptrandom.FillWeight();
        sig.TaggE = taggerhit.PhotonEnergy;
        sig.TaggT = taggerhit.Time;
        sig.TaggCh = taggerhit.Channel;

        // find best combination for each Tagger hit
        best_prob_fit = -std_ext::inf;
        best_comb_fit = cands.size();

        for (size_t i = 0; i < cands.size(); i++) {  // loop to test all different combinations
            // ensure the possible proton candidate is kinematically allowed
            if (std_ext::radian_to_degree(comb.back()->Theta) > 90.) {
                shift_right(comb);
                continue;
            }
            h.steps->Fill("proton #vartheta", 1);

            // require 2 PID entries for the eta candidate
            if (std::count_if(comb.begin(), comb.end()-1, [](TCandidatePtr c){ return c->VetoEnergy; }) < 2) {
                shift_right(comb);
                continue;
            }
            h.steps->Fill("2 PIDs", 1);

            photons.clear();
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
            etap.SetXYZT(0,0,0,0);
            for (size_t j = 0; j < comb.size()-1; j++) {
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, comb.at(j)));
                etap += TParticle(ParticleTypeDatabase::Photon, comb.at(j));
            }
            h.etapIM->Fill(etap.M(), sig.TaggW);

            // do the fitting and check if the combination is better than the previous best
            if (!doFit_checkProb(taggerhit, proton, photons, h, sig, best_prob_fit)) {
                shift_right(comb);
                continue;
            }

            best_comb_fit = i;

            shift_right(comb);
        }

        // only fill tree if a valid combination for the current Tagger hit was found
        if (best_comb_fit >= cands.size() || !isfinite(best_prob_fit))
            continue;

        sig.Tree->Fill();
        h.steps->Fill("Tree filled", 1);
    }

    if (best_comb_fit >= cands.size() || !isfinite(best_prob_fit))
        return;
    h.steps->Fill("best comb", 1);

    // restore combinations with best chi2
    //for (size_t i = 0; i < best_comb_fit; i++)
    while (best_comb_fit-- > 0)
        shift_right(comb);

    // sort the eta final state according to their Veto energies
    sort(comb.begin(), comb.end()-1,
         [] (const TCandidatePtr& a, const TCandidatePtr& b) {
            return a->VetoEnergy > b->VetoEnergy;
         });

    // do an anti pi0 cut on the combinations e+g and e-g
    // (assuming the photon deposited the least energy in the PIDs)
    if (ANTI_PI0_CUT) {
        const interval<double> pion_cut(ANTI_PI0_LOW, ANTI_PI0_HIGH);
        TLorentzVector pi0;
        const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};
        for (const auto pi0_comb : pi0_combs) {
            pi0 = TLorentzVector(0., 0., 0., 0.);
            for (const auto idx : pi0_comb)
                pi0 += TParticle(ParticleTypeDatabase::Photon, comb.at(idx));
            // apply an anti pion cut
            if (pion_cut.Contains(pi0.M()))
                return;
        }
        h.steps->Fill("anti #pi^{0} cut", 1);
    }

    proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
    etap.SetXYZT(0,0,0,0);
    for (size_t i = 0; i < comb.size()-1; i++)
        etap += TParticle(ParticleTypeDatabase::Photon, comb.at(i));
    h.etapIM_cand->Fill(etap.M());
    h_protonVeto->Fill(comb.back()->VetoEnergy);
    h_pTheta->Fill(std_ext::radian_to_degree(comb.back()->Theta));

    // at this point a possible eta Dalitz candidate was found, work only with eta final state
    comb.pop_back();

    const TCandidatePtr& l1 = comb.at(0);
    const TCandidatePtr& l2 = comb.at(1);
    // suppress conversion decays
    if (l1->FindVetoCluster()->CentralElement == l2->FindVetoCluster()->CentralElement)
        return;
    h.steps->Fill("distinct PID", 1);
    const double eeIM = (TParticle(ParticleTypeDatabase::eMinus, l1) + TParticle(ParticleTypeDatabase::eMinus, l2)).M();
    h_IM2d->Fill(etap.M(), eeIM);
    h.IM2d->Fill(etap.M(), eeIM);

    // test effective cluster radius to distinguish between leptons and charged pions
    double effective_radius = calc_effective_radius(l1);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l1->FindCaloCluster()->Energy, effective_radius);
    }
    effective_radius = calc_effective_radius(l2);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l2->FindCaloCluster()->Energy, effective_radius);
    }
    double lateral_moment = lat_moment(l1);
    if (isfinite(lateral_moment)) {
        h.lat_moment->Fill(lateral_moment);
        h.lat_moment_E->Fill(l1->FindCaloCluster()->Energy, lateral_moment);
    }
    lateral_moment = lat_moment(l2);
    if (isfinite(lateral_moment)) {
        h.lat_moment->Fill(lateral_moment);
        h.lat_moment_E->Fill(l2->FindCaloCluster()->Energy, lateral_moment);
    }

    // test cluster size compared to energy
    h.cluster_size->Fill(l1->ClusterSize);
    h.cluster_size->Fill(l2->ClusterSize);
    h.cluster_size_E->Fill(l1->FindCaloCluster()->Energy, l1->ClusterSize);
    h.cluster_size_E->Fill(l2->FindCaloCluster()->Energy, l2->ClusterSize);

    h.etapIM_final->Fill(etap.M());
    h_etapIM_final->Fill(etap.M());
    h.hCopl_final->Fill(std_ext::radian_to_degree(abs(etap.Phi() - proton->Phi())) - 180.);

    h_counts->Fill(decaystring.c_str(), 1);
}

void EtapDalitz::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << h_IM2d << endc;

    for (auto& entry : channels)
        entry.second.Show();

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
                                 const TParticlePtr proton,
                                 const TParticleList photons,
                                 PerChannel_t& h,
                                 SigTree_t& t,
                                 double& best_prob_fit)
{
    TLorentzVector etap(0,0,0,0);
    TLorentzVector etap_kinfit(0,0,0,0);
    TLorentzVector etap_kinfit_freeZ(0,0,0,0);
    TLorentzVector etap_treefit(0,0,0,0);
    TLorentzVector etap_treefit_freeZ(0,0,0,0);

    for (const auto& g : photons)
        etap += *g;

    /* kinematical checks to reduce computing time */
    const interval<double> coplanarity({-25, 25});
    const interval<double> mm = ParticleTypeDatabase::Proton.GetWindow(300);

    const double copl = std_ext::radian_to_degree(abs(etap.Phi() - proton->Phi())) - 180.;
    h.hCopl->Fill(copl, t.TaggW);
    if (!coplanarity.Contains(copl))
        return false;
    h.steps->Fill("coplanarity", 1);

    LorentzVec missing = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
    missing -= etap;
    h.MM->Fill(missing.M(), t.TaggW);
    if (!mm.Contains(missing.M()))
        return false;
    h.steps->Fill("missing mass", 1);


    /* now start with the kinematic fitting */

    // treefit
    APLCON::Result_t treefit_result;

    treefitter_etap.PrepareFits(taggerhit.PhotonEnergy, proton, photons);

    // works this way because only one combination needs to be fitted
    treefitter_etap.NextFit(treefit_result);

    if (USE_TREEFIT) {
        if (treefit_result.Status != APLCON::Result_Status_t::Success)
            return false;
        h.steps->Fill("treefit", 1);
    }

    // treefit free Z vertex
    APLCON::Result_t treefit_freeZ_result;

    treefitter_etap_freeZ.PrepareFits(taggerhit.PhotonEnergy, proton, photons);

    treefitter_etap_freeZ.NextFit(treefit_freeZ_result);


    // kinfit
    auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, proton, photons);

    if (!USE_TREEFIT) {
        if (kinfit_result.Status != APLCON::Result_Status_t::Success)
            return false;
        h.steps->Fill("kinfit", 1);
    }

    // kinfit free Z vertex
    auto kinfit_freeZ_result = kinfit_freeZ.DoFit(taggerhit.PhotonEnergy, proton, photons);


    const double treefit_chi2 = treefit_result.ChiSquare;
    const double treefit_prob = treefit_result.Probability;
    const int treefit_iterations = treefit_result.NIterations;
    const double treefit_freeZ_chi2 = treefit_freeZ_result.ChiSquare;
    const double treefit_freeZ_prob = treefit_freeZ_result.Probability;
    const int treefit_freeZ_iterations = treefit_freeZ_result.NIterations;
    const double kinfit_chi2 = kinfit_result.ChiSquare;
    const double kinfit_prob = kinfit_result.Probability;
    const int kinfit_iterations = kinfit_freeZ_result.NIterations;
    const double kinfit_freeZ_chi2 = kinfit_freeZ_result.ChiSquare;
    const double kinfit_freeZ_prob = kinfit_freeZ_result.Probability;
    const int kinfit_freeZ_iterations = kinfit_freeZ_result.NIterations;

    h.treefitChi2->Fill(treefit_chi2);
    h.treefitProb->Fill(treefit_prob);
    h.treefitIter->Fill(treefit_iterations);
    h.treefit_freeZ_chi2->Fill(treefit_freeZ_chi2);
    h.treefit_freeZ_prob->Fill(treefit_freeZ_prob);
    h.treefit_freeZ_iter->Fill(treefit_freeZ_iterations);
    h.kinfitChi2->Fill(kinfit_chi2);
    h.kinfitProb->Fill(kinfit_prob);
    h.kinfitIter->Fill(kinfit_iterations);
    h.kinfit_freeZ_chi2->Fill(kinfit_freeZ_chi2);
    h.kinfit_freeZ_prob->Fill(kinfit_freeZ_prob);
    h.kinfit_freeZ_iter->Fill(kinfit_freeZ_iterations);

    // determine which probability should be used to find the best candidate combination
    double prob = USE_TREEFIT ? treefit_prob : kinfit_prob;

    if (PROBABILITY_CUT) {
        if (prob < PROBABILITY)
            return false;
        h.steps->Fill("probability", 1);
    }

    h_etap->Fill(etap.E() - etap.M(), std_ext::radian_to_degree(etap.Theta()), t.TaggW);
    h_proton->Fill(proton->E - proton->M(), std_ext::radian_to_degree(proton->Theta()), t.TaggW);

    // check if a better probability has been found
    if (!std_ext::copy_if_greater(best_prob_fit, prob))
        return false;


    /* Gather relevant fitter information and update branches */

    // update branches with fitter and general particle information
    t.kinfit_chi2 = kinfit_chi2;
    t.kinfit_probability = kinfit_prob;
    t.kinfit_iterations = kinfit_iterations;
    t.kinfit_DoF = kinfit_result.NDoF;
    t.kinfit_freeZ_chi2 = kinfit_freeZ_chi2;
    t.kinfit_freeZ_probability = kinfit_freeZ_prob;
    t.kinfit_freeZ_iterations = kinfit_freeZ_iterations;
    t.kinfit_freeZ_DoF = kinfit_freeZ_result.NDoF;
    t.treefit_chi2 = treefit_chi2;
    t.treefit_probability = treefit_prob;
    t.treefit_iterations = treefit_iterations;
    t.treefit_DoF = treefit_result.NDoF;
    t.treefit_freeZ_chi2 = treefit_freeZ_chi2;
    t.treefit_freeZ_probability = treefit_freeZ_prob;
    t.treefit_freeZ_iterations = treefit_freeZ_iterations;
    t.treefit_freeZ_DoF = treefit_freeZ_result.NDoF;

    t.p               = *proton;
    t.p_Time          = proton->Candidate->Time;
    t.p_PSA           = getPSAVector(proton);
    t.p_vetoE         = proton->Candidate->VetoEnergy;
    t.p_detector      = getDetectorAsInt(proton->Candidate->Detector);
    t.p_clusterSize   = proton->Candidate->ClusterSize;
    t.p_centralElem   = proton->Candidate->FindCaloCluster()->CentralElement;
    t.p_effect_radius = calc_effective_radius(proton->Candidate);
    t.p_lat_moment    = lat_moment(proton->Candidate);
    t.p_vetoChannel   = -1;
    if (proton->Candidate->VetoEnergy)
        t.p_vetoChannel = proton->Candidate->FindVetoCluster()->CentralElement;

    for (size_t i = 0; i < N_FINAL_STATE-1; ++i) {
        t.photons().at(i)               = *(photons.at(i));
        t.photons_Time().at(i)          = photons.at(i)->Candidate->Time;
        t.photons_vetoE().at(i)         = photons.at(i)->Candidate->VetoEnergy;
        t.photons_PSA().at(i)           = getPSAVector(photons.at(i));
        t.photons_detector().at(i)      = getDetectorAsInt(photons.at(i)->Candidate->Detector);
        t.photons_clusterSize().at(i)   = photons.at(i)->Candidate->ClusterSize;
        t.photons_centralElem().at(i)   = photons.at(i)->Candidate->FindCaloCluster()->CentralElement;
        t.photons_effect_radius().at(i) = calc_effective_radius(photons.at(i)->Candidate);
        t.photons_lat_moment().at(i)    = lat_moment(photons.at(i)->Candidate);
        t.photons_vetoChannel().at(i)   = -1;
        if (photons.at(i)->Candidate->VetoEnergy)
            t.photons_vetoChannel().at(i) = photons.at(i)->Candidate->FindVetoCluster()->CentralElement;
    }

    t.etap = etap;
    t.mm = missing;
    t.copl = copl;

    // now handle the different fitted particle information separately

    // kinfit
    if (kinfit_result.Status == APLCON::Result_Status_t::Success) {
        auto kinfitted_photons = kinfit.GetFittedPhotons();
        auto kinfitted_proton = kinfit.GetFittedProton();
        auto kinfitted_beam = kinfit.GetFittedBeamE();
        auto kinfitted_beam_pull = kinfit.GetBeamEPull();
        auto kinfit_particles = kinfit.GetFitParticles();

        assert(kinfit_particles.size() == N_FINAL_STATE);

        for (const auto& g : kinfitted_photons)
            etap_kinfit += *g;
        h.etapIM_kinfit->Fill(etap_kinfit.M(), t.TaggW);

        h.kinfit_ZVertex->Fill(kinfit.GetFittedZVertex());

        // update tree branches
        t.beam_E_kinfitted = kinfitted_beam;
        t.beam_kinfit_E_pull = kinfitted_beam_pull;
        t.kinfit_ZVertex = kinfit.GetFittedZVertex();
        t.kinfit_ZVertex_pull = kinfit.GetZVertexPull();

        t.p_kinfitted     = *kinfitted_proton;

        t.p_kinfit_theta_pull        = kinfit_particles.at(0).GetPulls().at(1);
        t.p_kinfit_phi_pull          = kinfit_particles.at(0).GetPulls().at(2);

        for (size_t i = 0; i < N_FINAL_STATE-1; ++i) {
            t.photons_kinfitted().at(i)     = *(kinfitted_photons.at(i));

            t.photon_kinfit_E_pulls().at(i)            = kinfit_particles.at(i+1).GetPulls().at(0);
            t.photon_kinfit_theta_pulls().at(i)        = kinfit_particles.at(i+1).GetPulls().at(1);
            t.photon_kinfit_phi_pulls().at(i)          = kinfit_particles.at(i+1).GetPulls().at(2);
        }

        t.etap_kinfit = etap_kinfit;
    }

    // kinfit with free Z vertex
    if (kinfit_freeZ_result.Status == APLCON::Result_Status_t::Success) {
        auto kinfit_freeZ_photons = kinfit_freeZ.GetFittedPhotons();
        auto kinfit_freeZ_proton = kinfit_freeZ.GetFittedProton();
        auto kinfit_freeZ_beam = kinfit_freeZ.GetFittedBeamE();
        auto kinfit_freeZ_beam_pull = kinfit_freeZ.GetBeamEPull();
        auto kinfit_freeZ_particles = kinfit_freeZ.GetFitParticles();

        assert(kinfit_freeZ_particles.size() == N_FINAL_STATE);

        for (const auto& g : kinfit_freeZ_photons)
            etap_kinfit_freeZ += *g;
        h.etapIM_kinfit_freeZ->Fill(etap_kinfit_freeZ.M(), t.TaggW);

        h.kinfit_freeZ_ZVertex->Fill(kinfit_freeZ.GetFittedZVertex());

        // update tree branches
        t.beam_E_kinfit_freeZ = kinfit_freeZ_beam;
        t.beam_kinfit_freeZ_E_pull = kinfit_freeZ_beam_pull;
        t.kinfit_freeZ_ZVertex = kinfit_freeZ.GetFittedZVertex();
        t.kinfit_freeZ_ZVertex_pull = kinfit_freeZ.GetZVertexPull();

        t.p_kinfit_freeZ  = *kinfit_freeZ_proton;

        t.p_kinfit_freeZ_theta_pull  = kinfit_freeZ_particles.at(0).GetPulls().at(1);
        t.p_kinfit_freeZ_phi_pull    = kinfit_freeZ_particles.at(0).GetPulls().at(2);

        for (size_t i = 0; i < N_FINAL_STATE-1; ++i) {
            t.photons_kinfit_freeZ().at(i)  = *(kinfit_freeZ_photons.at(i));

            t.photon_kinfit_freeZ_E_pulls().at(i)      = kinfit_freeZ_particles.at(i+1).GetPulls().at(0);
            t.photon_kinfit_freeZ_theta_pulls().at(i)  = kinfit_freeZ_particles.at(i+1).GetPulls().at(1);
            t.photon_kinfit_freeZ_phi_pulls().at(i)    = kinfit_freeZ_particles.at(i+1).GetPulls().at(2);
        }

        t.etap_kinfit_freeZ = etap_kinfit_freeZ;
    }

    // treefit
    if (treefit_result.Status == APLCON::Result_Status_t::Success) {
        auto treefitted_photons = treefitter_etap.GetFittedPhotons();
        auto treefitted_proton = treefitter_etap.GetFittedProton();
        auto treefitted_beam = treefitter_etap.GetFittedBeamE();
        auto treefitted_beam_pull = treefitter_etap.GetBeamEPull();
        auto treefit_particles = treefitter_etap.GetFitParticles();

        assert(treefit_particles.size() == N_FINAL_STATE);

        for (const auto& g : treefitted_photons)
            etap_treefit += *g;
        h.etapIM_treefit->Fill(etap_treefit.M(), t.TaggW);

        h.treefit_ZVertex->Fill(treefitter_etap.GetFittedZVertex());

        // update tree branches
        t.beam_E_treefitted = treefitted_beam;
        t.beam_treefit_E_pull = treefitted_beam_pull;
        t.treefit_ZVertex = treefitter_etap.GetFittedZVertex();
        t.treefit_ZVertex_pull = treefitter_etap.GetZVertexPull();

        t.p_treefitted    = *treefitted_proton;

        t.p_treefit_theta_pull       = treefit_particles.at(0).GetPulls().at(1);
        t.p_treefit_phi_pull         = treefit_particles.at(0).GetPulls().at(2);

        for (size_t i = 0; i < N_FINAL_STATE-1; ++i) {
            t.photons_treefitted().at(i)    = *(treefitted_photons.at(i));

            t.photon_treefit_E_pulls().at(i)           = treefit_particles.at(i+1).GetPulls().at(0);
            t.photon_treefit_theta_pulls().at(i)       = treefit_particles.at(i+1).GetPulls().at(1);
            t.photon_treefit_phi_pulls().at(i)         = treefit_particles.at(i+1).GetPulls().at(2);
        }

        t.etap_treefit = etap_treefit;
    }

    // treefit with free Z vertex
    if (treefit_freeZ_result.Status == APLCON::Result_Status_t::Success) {
        auto treefit_freeZ_photons = treefitter_etap_freeZ.GetFittedPhotons();
        auto treefit_freeZ_proton = treefitter_etap_freeZ.GetFittedProton();
        auto treefit_freeZ_beam = treefitter_etap_freeZ.GetFittedBeamE();
        auto treefit_freeZ_beam_pull = treefitter_etap_freeZ.GetBeamEPull();
        auto treefit_freeZ_particles = treefitter_etap_freeZ.GetFitParticles();

        assert(treefit_freeZ_particles.size() == N_FINAL_STATE);

        for (const auto& g : treefit_freeZ_photons)
            etap_treefit_freeZ += *g;
        h.etapIM_treefit_freeZ->Fill(etap_treefit_freeZ.M(), t.TaggW);

        h.treefit_freeZ_ZVertex->Fill(treefitter_etap_freeZ.GetFittedZVertex());

        // update tree branches
        t.beam_E_treefit_freeZ = treefit_freeZ_beam;
        t.beam_treefit_freeZ_E_pull = treefit_freeZ_beam_pull;
        t.treefit_freeZ_ZVertex = treefitter_etap_freeZ.GetFittedZVertex();
        t.treefit_freeZ_ZVertex_pull = treefitter_etap_freeZ.GetZVertexPull();

        t.p_treefit_freeZ = *treefit_freeZ_proton;

        t.p_treefit_freeZ_theta_pull = treefit_freeZ_particles.at(0).GetPulls().at(1);
        t.p_treefit_freeZ_phi_pull   = treefit_freeZ_particles.at(0).GetPulls().at(2);

        for (size_t i = 0; i < N_FINAL_STATE-1; ++i) {
            t.photons_treefit_freeZ().at(i) = *(treefit_freeZ_photons.at(i));

            t.photon_treefit_freeZ_E_pulls().at(i)     = treefit_freeZ_particles.at(i+1).GetPulls().at(0);
            t.photon_treefit_freeZ_theta_pulls().at(i) = treefit_freeZ_particles.at(i+1).GetPulls().at(1);
            t.photon_treefit_freeZ_phi_pulls().at(i)   = treefit_freeZ_particles.at(i+1).GetPulls().at(2);
        }

        t.etap_treefit_freeZ = etap_treefit_freeZ;
    }

    return true;
}

void EtapDalitz::count_clusters(const TCandidateList& cands)
{
    size_t nCB = 0, nTAPS = 0;
    for (auto p : cands.get_iter())
        if (p->Detector & Detector_t::Type_t::CB)
            nCB++;
        else if (p->Detector & Detector_t::Type_t::TAPS)
            nTAPS++;
    h_cluster_CB->Fill(nCB);
    h_cluster_TAPS->Fill(nTAPS);
}

bool EtapDalitz::q2_preselection(const TEventData& data, const double threshold = 50.) const
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

EtapDalitz::ReactionChannel_t::~ReactionChannel_t()
{}

EtapDalitz::ReactionChannel_t::ReactionChannel_t(const string &n):
    name(n)
{}

EtapDalitz::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<EtapDalitz::decaytree_t> &t, const int c):
    name(utils::ParticleTools::GetDecayString(t)),
    tree(t),
    color(c)
{}

EtapDalitz::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<EtapDalitz::decaytree_t> &t, const string &n, const int c):
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
    const bool MC = event.Reconstructed().ID.isSet(TID::Flags_t::MC);

    triggersimu.ProcessEvent(event);
    const auto& cands = event.Reconstructed().Candidates;

    if (t->nCands != N_FINAL_STATE)
        return;

    // set fitter uncertainty models
    kinfit.SetUncertaintyModel(MC ? model_MC : model_data);
    treefitter_etap.SetUncertaintyModel(MC ? model_MC : model_data);

    TParticlePtr proton;
    TParticleList photons;

    if (!USE_KINFIT) {
        if (!simple2CB1TAPS(cands, proton, photons))
            return;

        for (const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {  // loop over all tagger hits
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

            while (treefitter_etap.NextFit(treefit_result))
                if (treefit_result.Status != APLCON::Result_Status_t::Success)
                    continue;

            // kinfit

            auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, proton, photons);

            // fill the tree with the fitted values
            fill_tree(treefit_result, kinfit_result, proton, photons);
            t->Tree->Fill();
        }

        return;
    }


    /* use kinematic fitting to determine the proton */
    TCandidatePtrList comb;
    for (auto p : cands.get_iter())
        comb.emplace_back(p);

    double best_prob_fit = -std_ext::inf;
    size_t best_comb_fit = cands.size();
    for (const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {  // loop over all tagger hits
        promptrandom->SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if (promptrandom->State() == PromptRandom::Case::Outside)
            continue;

        t->TaggW = promptrandom->FillWeight();
        t->TaggE = taggerhit.PhotonEnergy;
        t->TaggT = taggerhit.Time;
        t->TaggCh = taggerhit.Channel;

        // find best combination for each Tagger hit
        best_prob_fit = -std_ext::inf;
        best_comb_fit = cands.size();

        for (size_t i = 0; i < cands.size(); i++) {  // loop to test all different combinations
            // ensure the possible proton candidate is kinematically allowed
            if (std_ext::radian_to_degree(comb.back()->Theta) > 90.) {
                shift_right(comb);
                continue;
            }

            photons.clear();
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
            for (size_t j = 0; j < comb.size()-1; j++)
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, comb.at(j)));

            // do the fitting and check if the combination is better than the previous best
            if (!doFit_checkProb(taggerhit, proton, photons, best_prob_fit)) {
                shift_right(comb);
                continue;
            }

            best_comb_fit = i;

            shift_right(comb);
        }

        // only fill tree if a valid combination for the current Tagger hit was found
        if (best_comb_fit >= cands.size() || !isfinite(best_prob_fit))
            continue;

        t->Tree->Fill();
    }
}

void Etap2g::fill_tree(const APLCON::Result_t& treefit_result,
                       const APLCON::Result_t& kinfit_result,
                       const TParticlePtr proton,
                       const TParticleList photons)
{
    TLorentzVector etap;
    TLorentzVector etap_kinfit;
    TLorentzVector etap_treefit;

    auto sumlv = [] (TLorentzVector sum, TParticlePtr p) {
        return sum += *p;
    };

    auto kinfit_photons = kinfit.GetFittedPhotons();
    auto treefit_photons = treefitter_etap.GetFittedPhotons();
    auto kinfit_particles = kinfit.GetFitParticles();
    auto treefit_particles = treefitter_etap.GetFitParticles();

    etap = std::accumulate(photons.begin(), photons.end(), LorentzVec(), sumlv);
    etap_kinfit = std::accumulate(kinfit_photons.begin(), kinfit_photons.end(), LorentzVec(), sumlv);
    etap_treefit = std::accumulate(treefit_photons.begin(), treefit_photons.end(), LorentzVec(), sumlv);

    t->kinfit_chi2 = kinfit_result.ChiSquare;
    t->kinfit_probability = kinfit_result.Probability;
    t->kinfit_iterations = kinfit_result.NIterations;
    t->kinfit_DoF = kinfit_result.NDoF;
    t->treefit_chi2 = treefit_result.ChiSquare;
    t->treefit_probability = treefit_result.Probability;
    t->treefit_iterations = treefit_result.NIterations;
    t->treefit_DoF = treefit_result.NDoF;

    t->beam_E_kinfitted = kinfit.GetFittedBeamE();
    t->beam_kinfit_E_pull = kinfit.GetBeamEPull();
    t->beam_E_treefitted = treefitter_etap.GetFittedBeamE();
    t->beam_treefit_E_pull = treefitter_etap.GetBeamEPull();
    t->kinfit_ZVertex = kinfit.GetFittedZVertex();
    t->kinfit_ZVertex_pull = kinfit.GetZVertexPull();
    t->treefit_ZVertex = treefitter_etap.GetFittedZVertex();
    t->treefit_ZVertex_pull = treefitter_etap.GetZVertexPull();

    t->p               = *proton;
    t->p_kinfitted     = *(kinfit.GetFittedProton());
    t->p_treefitted    = *(treefitter_etap.GetFittedProton());
    t->p_Time          = proton->Candidate->Time;
    t->p_PSA           = getPSAVector(proton);
    t->p_vetoE         = proton->Candidate->VetoEnergy;
    t->p_detector      = getDetectorAsInt(proton->Candidate->Detector);
    t->p_clusterSize   = proton->Candidate->ClusterSize;
    t->p_centralElem   = proton->Candidate->FindCaloCluster()->CentralElement;
    t->p_vetoChannel   = -1;
    if (proton->Candidate->VetoEnergy)
        t->p_vetoChannel = proton->Candidate->FindVetoCluster()->CentralElement;

    t->p_kinfit_theta_pull        = kinfit_particles.at(0).GetPulls().at(1);
    t->p_kinfit_phi_pull          = kinfit_particles.at(0).GetPulls().at(2);
    t->p_treefit_theta_pull       = treefit_particles.at(0).GetPulls().at(1);
    t->p_treefit_phi_pull         = treefit_particles.at(0).GetPulls().at(2);

    for (size_t i = 0; i < N_FINAL_STATE-1; ++i) {
        t->photons().at(i)               = *(photons.at(i));
        t->photons_kinfitted().at(i)     = *(kinfit_photons.at(i));
        t->photons_treefitted().at(i)    = *(treefit_photons.at(i));
        t->photons_Time().at(i)          = photons.at(i)->Candidate->Time;
        t->photons_vetoE().at(i)         = photons.at(i)->Candidate->VetoEnergy;
        t->photons_detector().at(i)      = getDetectorAsInt(photons.at(i)->Candidate->Detector);
        t->photons_centralElem().at(i)   = photons.at(i)->Candidate->FindCaloCluster()->CentralElement;
        t->photons_vetoChannel().at(i)   = -1;
        if (photons.at(i)->Candidate->VetoEnergy)
            t->photons_vetoChannel().at(i) = photons.at(i)->Candidate->FindVetoCluster()->CentralElement;

        t->photon_kinfit_E_pulls().at(i)            = kinfit_particles.at(i+1).GetPulls().at(0);
        t->photon_kinfit_theta_pulls().at(i)        = kinfit_particles.at(i+1).GetPulls().at(1);
        t->photon_kinfit_phi_pulls().at(i)          = kinfit_particles.at(i+1).GetPulls().at(2);
        t->photon_treefit_E_pulls().at(i)           = treefit_particles.at(i+1).GetPulls().at(0);
        t->photon_treefit_theta_pulls().at(i)       = treefit_particles.at(i+1).GetPulls().at(1);
        t->photon_treefit_phi_pulls().at(i)         = treefit_particles.at(i+1).GetPulls().at(2);
    }

    t->etap = etap;
    t->etap_kinfit = etap_kinfit;
    t->etap_treefit = etap_treefit;
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
                             const TParticleList photons,
                             double& best_prob_fit)
{
    TLorentzVector etap(0,0,0,0);

    for (const auto& g : photons)
        etap += *g;

    /* kinematical checks to reduce computing time */
    const interval<double> coplanarity({-25, 25});
    const interval<double> mm = ParticleTypeDatabase::Proton.GetWindow(300);

    const double copl = std_ext::radian_to_degree(abs(etap.Phi() - proton->Phi())) - 180.;
    if (!coplanarity.Contains(copl))
        return false;

    LorentzVec missing = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
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
    double prob = USE_TREEFIT ? treefit_result.Probability : kinfit_result.Probability;

//    if (PROBABILITY_CUT) {
//        if (prob < PROBABILITY)
//            return false;
//        h.steps->Fill("probability", 1);
//    }

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
