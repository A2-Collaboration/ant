#include "PID_Energy_etaDalitz.h"

#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

template<typename T>
void PID_Energy_etaDalitz::shift_right(std::vector<T>& v)
{
    std::rotate(v.begin(), v.end() -1, v.end());
}

void PID_Energy_etaDalitz::remove_char(std::string& str, char ch)
{
    str.erase(std::remove(str.begin(), str.end(), ch), str.end());
}

void PID_Energy_etaDalitz::remove_chars(std::string& str, std::initializer_list<char> chars)
{
    for (const auto ch : chars)
        remove_char(str, ch);
}

double PID_Energy_etaDalitz::calc_effective_radius(const TCandidatePtr cand)
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

ParticleTypeTree PID_Energy_etaDalitz::base_tree()
{
    ParticleTypeTree t = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(ParticleTypeDatabase::BeamProton);
    t->CreateDaughter(ParticleTypeDatabase::Proton);
    return t;
}

ParticleTypeTree PID_Energy_etaDalitz::eta_3g()
{
    auto t = base_tree();
    auto eta = t->CreateDaughter(ParticleTypeDatabase::Eta);
    eta->CreateDaughter(ParticleTypeDatabase::Photon);
    eta->CreateDaughter(ParticleTypeDatabase::Photon);
    eta->CreateDaughter(ParticleTypeDatabase::Photon);
    return t;
}

APLCON::Fit_Settings_t PID_Energy_etaDalitz::MakeFitSettings(unsigned max_iterations)
{
    auto settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = max_iterations;
    return settings;
}

PID_Energy_etaDalitz::Tree_t::Tree_t()
{}

PID_Energy_etaDalitz::PerChannel_t::PerChannel_t(const std::string& Name, const string& Title, HistogramFactory& hf):
    title(Title),
    name(Name)
{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);
    const BinSettings pid_channels(detector->GetNChannels());
    const BinSettings veto_energy(1000, 0, 10);
    const BinSettings energy(1200);

    eegPID = hf.makeTH2D(title + " PID 2 charged 1 neutral", "PID Energy [MeV]", "#", veto_energy, pid_channels, name + "_eegPID");
    steps = hf.makeTH1D(title + " Accepted Events", "step", "#", BinSettings(10), name + "_steps");
    etaIM = hf.makeTH1D(title + " IM #eta all comb", "IM [MeV]", "#", energy, name + "_etaIM");
    etaIM_fit = hf.makeTH1D(title + " IM #eta fitted", "IM [MeV]", "#", energy, name + "_etaIM_fit");
    etaIM_cand = hf.makeTH1D(title + " IM #eta candidates", "IM [MeV]", "#", energy, name + "_etaIM_cand");
    etaIM_final = hf.makeTH1D(title + " IM #eta final", "IM [MeV]", "#", energy, name + "_etaIM_final");
    IM2d = hf.makeTH2D(title + " IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), name + "_IM2d");
    MM = hf.makeTH1D(title + " Missing Mass proton", "MM [MeV]", "#", BinSettings(1600), name + "_MM");
    hCopl = hf.makeTH1D(title + " Coplanarity #eta - proton all comb", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), name + "_hCopl");
    hCopl_final = hf.makeTH1D(title + " Coplanarity #eta - proton final", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), name + "_hCopl_final");
    hChi2 = hf.makeTH1D(title + " #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + "_hChi2");
    hProb = hf.makeTH1D(title + " Probability", "probability", "#", BinSettings(500, 0, 1), name + "_hProb");
    hIter = hf.makeTH1D(title + " # Iterations", "#iterations", "#", BinSettings(20), name + "_hIter");
    effect_rad = hf.makeTH1D(title + " Effective Radius", "R", "#", BinSettings(500, 0, 50), name + "effect_rad");
    effect_rad_E = hf.makeTH2D(title + " Effective Radius vs. Cluster Energy", "E [MeV]", "R", energy, BinSettings(500, 0, 50), name + "effect_rad_E");

    proton_E_theta = hf.makeTH2D(title + " proton", "E [MeV]", "#vartheta [#circ]", energy, BinSettings(360, 0, 180), name + "_e_theta");
}

void PID_Energy_etaDalitz::PerChannel_t::Show()
{
    //canvas("Per Channel: " + title) << drawoption("colz") << eegPID << endc;
    canvas("Per Channel: " + title) << steps
                                    << etaIM_fit
                                    << etaIM_final
                                    << hCopl_final
                                    << hChi2
                                    << hProb
                                    << endc;
}

void PID_Energy_etaDalitz::PerChannel_t::Fill(const TEventData& d)
{
    const auto& protons = d.Particles.Get(ParticleTypeDatabase::Proton);
    if (!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), p->Theta()*TMath::RadToDeg());
    }
}

PID_Energy_etaDalitz::PID_Energy_etaDalitz(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    kinfit("kinfit", 3,
           make_shared<uncertainty_model_t>(), false, MakeFitSettings(20)
           ),
    treefitter_eta("treefitter_eta", eta_3g(),
                   make_shared<uncertainty_model_t>(), false, {}, MakeFitSettings(20)
                   )
{
    promptrandom.AddPromptRange({-5, 5});
    promptrandom.AddRandomRange({-30, -10});
    promptrandom.AddRandomRange({10, 30});

    cb = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);

    t.CreateBranches(HistFac.makeTTree("tree"));

    const BinSettings pid_channels(detector->GetNChannels());
    const BinSettings energybins(1000, 0, 10);

    h_eegPID = HistFac.makeTH2D("PID 2 charged 1 neutral", "PID Energy [MeV]", "#",
                                energybins, pid_channels, "h_eegPID");
    h_counts = HistFac.makeTH1D("Events per Channel", "channel", "#", BinSettings(20), "h_counts");
    h_pTheta = HistFac.makeTH1D("#vartheta proton candidate", "#vartheta_{p} [#circ]", "#", BinSettings(720, 0, 180), "h_pTheta");
    h_protonVeto = HistFac.makeTH1D("Veto energy identified proton", "Veto [MeV]", "#", energybins, "h_protonVeto");
    h_etaIM_final = HistFac.makeTH1D("IM #eta final", "IM [MeV]", "#", BinSettings(1200), "h_etaIM_final");
    h_IM2d = HistFac.makeTH2D("IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), "h_IM2d");
    h_eta = HistFac.makeTH2D("Kinematics #eta", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(360, 0, 180), "h_eta");
    h_proton = HistFac.makeTH2D("Kinematics p", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(160, 0, 80), "h_proton");
}

void PID_Energy_etaDalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    const bool MC = event.HasMCTrue();

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

    //const auto& cands = event.Reconstructed().Candidates;
    const auto& data = event.Reconstructed();
    const auto& cands = data.Candidates;
    //const auto nCandidates = cands.size();
    t.nCands = cands.size();
    h.steps->Fill("seen", 1);

    if (cands.size() != 4)
        return;
    h.steps->Fill("#cands", 1);

    TLorentzVector eta;
    LorentzVec proton;
    //const interval<double> eta_im({ETA_IM-ETA_SIGMA, ETA_IM+ETA_SIGMA});
    const interval<double> coplanarity({-25, 25});
    TCandidatePtrList comb;
    for (auto p : cands.get_iter())
        comb.emplace_back(p);
/* old routine w/o kinfit
    bool found;
    for (size_t j = 0; j < cands.size(); j++) {  // loop to test all different combinations
        // ensure the possible proton candidate is kinematically allowed
        if (std_ext::radian_to_degree(comb.back()->Theta) > 60.)
            continue;
        found = true;
        eta.SetXYZT(0,0,0,0);
        for (size_t i = 0; i < comb.size()-1; i++)
            eta += TParticle(ParticleTypeDatabase::Photon, comb.at(i));
        proton = TParticle(ParticleTypeDatabase::Proton, comb.back());
        const double copl = std_ext::radian_to_degree(abs(eta.Phi() - proton.Phi())) - 180.;
        etaIM->Fill(eta.M());
        hCopl->Fill(copl);
        if (eta_im.Contains(eta.M()) && coplanarity.Contains(copl) && comb.back()->VetoEnergy
                && std::count_if(comb.begin(), comb.end()-1, [](TCandidatePtr c){ return c->VetoEnergy; }) >= 2)
            break;
        found = false;
        shift_right(comb);
    }
    if (!found)
        return;
*/

    // require at least 2 candidates with PID/Veto entries
    if (std::count_if(comb.begin(), comb.end(), [](TCandidatePtr c){ return c->VetoEnergy; }) < 2)
        return;
    h.steps->Fill("#Veto", 1);

    t.CBAvgTime = event.Reconstructed().Trigger.CBTiming;
    if (MC)
        t.CBAvgTime = 0;
    if (!isfinite(t.CBAvgTime))
        return;

    TParticleList photons;
    LorentzVec missing;
    std::vector<utils::Fitter::FitParticle> fitparticles;
    const interval<double> mm({ParticleTypeDatabase::Proton.Mass()-150., ParticleTypeDatabase::Proton.Mass()+150.});
    double best_prob = -std_ext::inf;
    size_t best_comb = cands.size();
    for (const TTaggerHit& taggerhit : data.TaggerHits) {  // loop over all tagger hits
        promptrandom.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h.steps->Fill("time window", 1);

        t.TaggW = promptrandom.FillWeight();
        t.TaggE = taggerhit.PhotonEnergy;
        t.TaggT = taggerhit.Time;
        t.TaggCh = taggerhit.Channel;

        for (size_t i = 0; i < cands.size(); i++) {  // loop to test all different combinations
            // ensure the possible proton candidate is kinematically allowed
            if (std_ext::radian_to_degree(comb.back()->Theta) > 60.) {
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
            proton = TParticle(ParticleTypeDatabase::Proton, comb.back());
            eta.SetXYZT(0,0,0,0);
            for (size_t j = 0; j < comb.size()-1; j++) {
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, comb.at(j)));
                eta += TParticle(ParticleTypeDatabase::Photon, comb.at(j));
            }
            t.eta = eta;
            h.etaIM->Fill(eta.M(), t.TaggW);

            const double copl = std_ext::radian_to_degree(abs(eta.Phi() - proton.Phi())) - 180.;
            t.copl = copl;
            h.hCopl->Fill(copl, t.TaggW);
            if (!coplanarity.Contains(copl)) {
                shift_right(comb);
                continue;
            }
            h.steps->Fill("coplanarity", 1);

            missing = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
            missing -= eta;
            t.missing_momentum = missing;
            h.MM->Fill(missing.M(), t.TaggW);
            if (!mm.Contains(missing.M())) {
                shift_right(comb);
                continue;
            }
            h.steps->Fill("missing mass", 1);

            APLCON::Result_t r;
            TParticleList fitted_photons;

            if (USE_TREEFIT) {
                treefitter_eta.SetEgammaBeam(taggerhit.PhotonEnergy);
                const auto& p = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
                treefitter_eta.SetProton(p);
                treefitter_eta.SetPhotons(photons);

                // works this way because only one combination needs to be fitted
                while (treefitter_eta.NextFit(r))
                    if(r.Status != APLCON::Result_Status_t::Success)
                        continue;

                if (r.Status != APLCON::Result_Status_t::Success) {
                    shift_right(comb);
                    continue;
                }
                h.steps->Fill("treefit", 1);

                fitted_photons = treefitter_eta.GetFittedPhotons();
            } else {
                kinfit.SetEgammaBeam(taggerhit.PhotonEnergy);
                const auto& p = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
                kinfit.SetProton(p);
                kinfit.SetPhotons(photons);

                r = kinfit.DoFit();

                if (r.Status != APLCON::Result_Status_t::Success) {
                    shift_right(comb);
                    continue;
                }
                h.steps->Fill("kinfit", 1);

                fitted_photons = kinfit.GetFittedPhotons();
            }

            const double chi2 = r.ChiSquare;
            const double prob = r.Probability;
            const int iterations = r.NIterations;
            h.hChi2->Fill(chi2);
            h.hProb->Fill(prob);
            h.hIter->Fill(iterations);

            if (prob < .05) {
                shift_right(comb);
                continue;
            }
            h.steps->Fill("probability", 1);

            h_eta->Fill(eta.E() - eta.M(), std_ext::radian_to_degree(eta.Theta()), t.TaggW);
            h_proton->Fill(proton.E - proton.M(), std_ext::radian_to_degree(proton.Theta()), t.TaggW);

            eta.SetXYZT(0,0,0,0);
            for (const auto& g : fitted_photons)
                eta += *g;
            t.eta_fit = eta;
            h.etaIM_fit->Fill(eta.M(), t.TaggW);

            if (prob > best_prob) {
                best_prob = prob;
                best_comb = i;
            }
            t.chi2 = chi2/r.NDoF;
            t.probability = prob;
            t.iterations = iterations;

            shift_right(comb);
        }
    }

    if (best_comb >= cands.size() || !isfinite(best_prob))
        return;
    h.steps->Fill("best comb", 1);

    // restore combinations with best chi2
    while (best_comb-- > 0)
        shift_right(comb);

    proton = TParticle(ParticleTypeDatabase::Proton, comb.back());
    eta.SetXYZT(0,0,0,0);
    for (size_t i = 0; i < comb.size()-1; i++)
        eta += TParticle(ParticleTypeDatabase::Photon, comb.at(i));
    h.etaIM_cand->Fill(eta.M());
    h_protonVeto->Fill(comb.back()->VetoEnergy);
    h_pTheta->Fill(std_ext::radian_to_degree(comb.back()->Theta));
    // at this point a possible eta Dalitz candidate was found, work only with eta final state
    comb.pop_back();

    sort(comb.begin(), comb.end(),
         [] (const TCandidatePtr& a, const TCandidatePtr& b) {
            return a->VetoEnergy > b->VetoEnergy;
         });

    const TCandidatePtr& l1 = comb.at(0);
    const TCandidatePtr& l2 = comb.at(1);
    // suppress conversion decays
    if (l1->FindVetoCluster()->CentralElement == l2->FindVetoCluster()->CentralElement)
        return;
    h.steps->Fill("distinct PID", 1);
    const double eeIM = (TParticle(ParticleTypeDatabase::eMinus, l1) + TParticle(ParticleTypeDatabase::eMinus, l2)).M();
    h_IM2d->Fill(eta.M(), eeIM);
    h.IM2d->Fill(eta.M(), eeIM);
    // suppress pi0
//    if (eeIM < 130.)
//        return;
//    h.steps->Fill("above #pi^{0}", 1);

    // test effective cluster radius to distinguish between leptons and charged pions
    double effective_radius = calc_effective_radius(l1);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l1->FindCaloCluster()->Energy, effective_radius);
    }
    effective_radius = calc_effective_radius(l2);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l1->FindCaloCluster()->Energy, effective_radius);
    }

    h.etaIM_final->Fill(eta.M());
    h_etaIM_final->Fill(eta.M());
    h.hCopl_final->Fill(std_ext::radian_to_degree(abs(eta.Phi() - proton.Phi())) - 180.);
    for (const TCandidatePtr& c : comb)
        if (c->VetoEnergy) {
            h.eegPID->Fill(c->VetoEnergy, c->FindVetoCluster()->CentralElement);
            h_eegPID->Fill(c->VetoEnergy, c->FindVetoCluster()->CentralElement);
            //h_eegPID->Fill(c->VetoEnergy*sin(c->Theta), c->FindVetoCluster()->CentralElement);
        }
    h_counts->Fill(decaystring.c_str(), 1);
    t.Tree->Fill();
}

void PID_Energy_etaDalitz::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << h_eegPID << endc;

    for (auto& entry : channels)
        entry.second.Show();
}

AUTO_REGISTER_PHYSICS(PID_Energy_etaDalitz)
