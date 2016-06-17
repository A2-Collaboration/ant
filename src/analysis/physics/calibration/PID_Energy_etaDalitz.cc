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

APLCON::Fit_Settings_t PID_Energy_etaDalitz::MakeFitSettings(unsigned max_iterations)
{
    auto settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = max_iterations;
    return settings;
}

PID_Energy_etaDalitz::Tree_t::Tree_t()
{}

PID_Energy_etaDalitz::PerChannel_t::PerChannel_t(const string& Title, HistogramFactory& hf):
    title(Title)
{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);
    const BinSettings pid_channels(detector->GetNChannels());
    const BinSettings veto_energy(1000, 0, 10);
    const BinSettings energy(1200);

    eegPID = hf.makeTH2D(title + " PID 2 charged 1 neutral", "PID Energy [MeV]", "#", veto_energy, pid_channels, title + "_eegPID");
    steps = hf.makeTH1D(title + " Accepted Events", "step", "#", BinSettings(10), title + "_steps");
    etaIM = hf.makeTH1D(title + " #eta IM all comb", "IM [MeV]", "#", energy, title + "_etaIM");
    etaIM_fit = hf.makeTH1D(title + " #eta IM fitted", "IM [MeV]", "#", energy, title + "_etaIM_fit");
    etaIM_cand = hf.makeTH1D(title + " #eta IM candidates", "IM [MeV]", "#", energy, title + "_etaIM_cand");
    etaIM_final = hf.makeTH1D(title + " #eta IM final", "IM [MeV]", "#", energy, title + "_etaIM_final");
    MM = hf.makeTH1D(title + " Missing Mass proton", "MM [MeV]", "#", BinSettings(1600), title + "_MM");
    hCopl = hf.makeTH1D(title + " Coplanarity #eta - proton all comb", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), title + "_hCopl");
    hCopl_final = hf.makeTH1D(title + " Coplanarity #eta - proton final", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), title + "_hCopl_final");
    hChi2 = hf.makeTH1D(title + " #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), title + "_hChi2");
    hProb = hf.makeTH1D(title + " Probability", "probability", "#", BinSettings(500, 0, 1), title + "_hProb");
    hIter = hf.makeTH1D(title + " # Iterations", "#iterations", "#", BinSettings(20), title + "_hIter");

    proton_E_theta = hf.makeTH2D(title, "E [MeV]", "#vartheta [#circ]", energy, BinSettings(360, 0, 180), title + "_e_theta");
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
           make_shared<uncertainty_model_t>(),
           PID_Energy_etaDalitz::MakeFitSettings(20)
           )
{
    promptrandom.AddPromptRange({-5, 5});
    promptrandom.AddRandomRange({-30, -10});
    promptrandom.AddRandomRange({10, 30});

    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);

    t.CreateBranches(HistFac.makeTTree("tree"));

    const BinSettings pid_channels(detector->GetNChannels());
    const BinSettings energybins(1000, 0, 10);

    h_eegPID = HistFac.makeTH2D("PID 2 charged 1 neutral", "PID Energy [MeV]", "#",
                                energybins, pid_channels, "h_eegPID");
    h_counts = HistFac.makeTH1D("Events per Channel", "channel", "#", BinSettings(20), "h_counts");
    h_pTheta = HistFac.makeTH1D("#vartheta proton candidate", "#vartheta_{p} [#circ]", "#", BinSettings(720, 0, 180), "h_pTheta");
    h_protonVeto = HistFac.makeTH1D("Veto energy identified proton", "Veto [MeV]", "#", energybins, "h_protonVeto");
    h_eta = HistFac.makeTH2D("Kinematics #eta", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(360, 0, 180), "h_eta");
    h_proton = HistFac.makeTH2D("Kinematics p", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(160, 0, 80), "h_proton");
}

void PID_Energy_etaDalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    const bool MC = event.HasMCTrue();

    std::string decaystring;
    if (MC)
        decaystring = utils::ParticleTools::GetProductionChannelString(event.MCTrue().ParticleTree);
    else
        decaystring = "data";

    auto c = channels.find(decaystring);
    if (c == channels.end())
        channels.insert({decaystring, PerChannel_t(decaystring, HistFac)});

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
    TLorentzVector proton;
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
    TLorentzVector missing;
    std::vector<utils::Fitter::FitParticle> fitparticles;
    const interval<double> mm({ParticleTypeDatabase::Proton.Mass()-150., ParticleTypeDatabase::Proton.Mass()+150.});
    double min_chi2 = std_ext::inf;
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

            kinfit.SetEgammaBeam(taggerhit.PhotonEnergy);
            const auto& p = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
            kinfit.SetProton(p);
            kinfit.SetPhotons(photons);

            auto result = kinfit.DoFit();

            if (result.Status != APLCON::Result_Status_t::Success) {
                shift_right(comb);
                continue;
            }
            h.steps->Fill("kinfit", 1);

            const double chi2 = result.ChiSquare;
            const double prob = result.Probability;
            const int iterations = result.NIterations;
            h.hChi2->Fill(chi2);
            h.hProb->Fill(prob);
            h.hIter->Fill(iterations);

            if (prob < .1) {
                shift_right(comb);
                continue;
            }
            h.steps->Fill("probability", 1);

            h_eta->Fill(eta.E() - eta.M(), std_ext::radian_to_degree(eta.Theta()), t.TaggW);
            h_proton->Fill(proton.E() - proton.M(), std_ext::radian_to_degree(proton.Theta()), t.TaggW);

            auto fitted_photons = kinfit.GetFittedPhotons();
            eta.SetXYZT(0,0,0,0);
            for (const auto& g : fitted_photons)
                eta += *g;
            t.eta_fit = eta;
            h.etaIM_fit->Fill(eta.M(), t.TaggW);

            if (chi2 < min_chi2) {
                min_chi2 = chi2;
                best_comb = i;
            }
            t.chi2 = chi2/result.NDoF;
            t.probability = prob;
            t.iterations = iterations;

            shift_right(comb);
        }
    }

    if (best_comb >= cands.size() || !isfinite(min_chi2))
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
    // suppress pi0
    if (eeIM > 130.)
        return;
    h.steps->Fill("below #pi^{0}", 1);

    h.etaIM_final->Fill(eta.M());
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
