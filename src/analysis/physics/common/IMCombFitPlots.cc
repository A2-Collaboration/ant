#include "IMCombFitPlots.h"
#include "base/Logger.h"
#include "analysis/utils/uncertainties/FitterSergey.h"

#include "TH1D.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;


template<typename T>
bool IMCombFitPlots::shift_right(std::vector<T>& v)
{
    rotate(v.begin(), v.end()-1, v.end());
    return true;
}

template<typename iter>
LorentzVec IMCombFitPlots::sumlv(iter start, iter end) {
    LorentzVec s;
    while (start != end) {
        s += **(start);
        ++start;
    }
    return s;
}

APLCON::Fit_Settings_t IMCombFitPlots::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    return settings;
}

const ParticleTypeTree IMCombFitPlots::ptreeSignal = ParticleTypeTreeDatabase::Get(signal);

IMCombFitPlots::IMCombFitPlots(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    MAX_GAMMA(opts->Get<unsigned>("MaxGamma", 5)),
    USE_MC_SIGNAL(opts->Get<bool>("Signal", false)),
    raw_2(MaxNGamma()-1, {prs}),
    raw_n(MaxNGamma()-1, {prs}),
    fit_2(MaxNGamma()-1, {prs}),
    fit_n(MaxNGamma()-1, {prs}),
    raw_mm(MaxNGamma()-1, {prs}),
    raw_pE(MaxNGamma()-1, {prs}),
    fit_pE(MaxNGamma()-1, {prs}),
    dEvE_all(MaxNGamma()-1, {prs}),
    dEvE(MaxNGamma()-1, {24, {prs}}),
    dEvE_all_combined({prs}),
    dEvE_combined(24, {prs}),
    model(make_shared<utils::UncertaintyModels::FitterSergey>())
{
    prs.AddPromptRange({-3,2});
    prs.AddRandomRange({-35,-10});
    prs.AddRandomRange({10,35});

    LOG(INFO) << "Promt Random Ratio = " << prs.Ratio();
    LOG(INFO) << "Test events with up to " << MaxNGamma() << " Photons";
    if (USE_MC_SIGNAL)
        LOG(INFO) << "For MC only the specified signal channel will be considered";

    unsigned n = MinNGamma()-1;
    while (++n <= MaxNGamma())
        kinfit.emplace_back(utils::KinFitter(model, opts->HasOption("SigmaZ"), MakeFitSettings(20)));

    if (opts->HasOption("SigmaZ")) {
        double sigma_z = 0.;
        std_ext::copy_if_greater(sigma_z, opts->Get<double>("SigmaZ", 0.));
        LOG(INFO) << "Fit Z vertex enabled with sigma = " << sigma_z;
        for (auto& fit : kinfit)
            fit.SetZVertexSigma(sigma_z);
    }

    auto pid_hists = HistogramFactory("PIDProtonHistograms", HistFac);

    dEvE_all_combined.MakeHistograms(pid_hists,"dEvE_all_combined","M+ combined all channels dEvE proton fitted",
                                           BinSettings(1200),BinSettings(100,0,10),"E_{p} [MeV]","E_{PID} [MeV]");
    for (size_t i = 0; i < MaxNGamma()-1; ++i) {
        raw_2.at(i).MakeHistograms(HistFac,"IM_raw_2_"+to_string(i+2),to_string(i+2)+" #gamma IM",BinSettings(1200),"2#gamma IM [MeV]","");
        raw_n.at(i).MakeHistograms(HistFac,"IM_raw_n_"+to_string(i+2),to_string(i+2)+" #gamma IM",BinSettings(1200),to_string(i+2)+"#gamma IM [MeV]","");
        fit_2.at(i).MakeHistograms(HistFac,"IM_fit_2_"+to_string(i+2),to_string(i+2)+" #gamma IM fitted",BinSettings(1200),"2#gamma IM [MeV]","");
        fit_n.at(i).MakeHistograms(HistFac,"IM_fit_n_"+to_string(i+2),to_string(i+2)+" #gamma IM fitted",BinSettings(1200),to_string(i+2)+"#gamma IM [MeV]","");
        raw_mm.at(i).MakeHistograms(pid_hists,"MM_raw_p_"+to_string(i+2),to_string(i+2)+" proton MM",BinSettings(1500),to_string(i+2)+"proton MM [MeV]","");
        raw_pE.at(i).MakeHistograms(pid_hists,"pE_raw_p_"+to_string(i+2),to_string(i+2)+" proton Energy",BinSettings(1200),to_string(i+2)+"E_{p} [MeV]","");
        fit_pE.at(i).MakeHistograms(pid_hists,"pE_fit_p_"+to_string(i+2),to_string(i+2)+" proton Energy fitted",BinSettings(1200),to_string(i+2)+"E_{p} [MeV]","");
        dEvE_all.at(i).MakeHistograms(pid_hists,"dEvE_all_"+to_string(i+2),to_string(i+2)+" all channels dEvE proton fitted",
                                      BinSettings(1200),BinSettings(100,0,10),"E_{p} [MeV]","E_{PID} [MeV]");
        unsigned j = 0;
        for (auto& hist : dEvE.at(i)) {
            hist.MakeHistograms(pid_hists,"dEvE_fit_p_"+to_string(i+2)+"_chan"+to_string(j),to_string(i+2)+" channel "+to_string(j)+" dEvE proton fitted",
                                BinSettings(1200),BinSettings(100,0,10),"E_{p} [MeV]","E_{PID} [MeV]");
            j++;
        }
    }
    unsigned i = 0;
    for (auto& hist : dEvE_combined) {
        hist.MakeHistograms(pid_hists,"dEvE_fit_p_combined_chan"+to_string(i),"M+ combined channel "+to_string(i)+" dEvE proton fitted",
                            BinSettings(1200),BinSettings(100,0,10),"E_{p} [MeV]","E_{PID} [MeV]");
        i++;
    }

    projections = pid_hists.makeTH2D("Projections of High Energy Tail of Protons", "E_{PID} [MeV]", "PID Channel",
                                   BinSettings(100,0,10), BinSettings(24), "projections_hep");
}

void IMCombFitPlots::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    const auto& cands = event.Reconstructed().Candidates;
    if (cands.size() > MaxNGamma()+1 || cands.size() < MinNGamma()+1)
        return;

    if (USE_MC_SIGNAL)
        if (event.Reconstructed().ID.isSet(TID::Flags_t::MC) && event.MCTrue().ParticleTree
            && !event.MCTrue().ParticleTree->IsEqual(ptreeSignal, utils::ParticleTools::MatchByParticleName))
            return;

    TCandidatePtrList comb;
    TParticleList fitted_photons;
    TParticlePtr fitted_proton;

    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {
        prs.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        comb.clear();
        for (auto p : cands.get_iter())
            comb.emplace_back(p);
        if (!find_best_comb(taggerhit, comb, fitted_photons, fitted_proton))
            continue;

        // proton and photons found, only photons needed
        TParticlePtr proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
        LorentzVec missing = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
        comb.pop_back();
        TParticleList photons;
        for (const auto& p : comb)
            photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, p));

        const LorentzVec sum_raw = sumlv(photons.begin(), photons.end());
        const LorentzVec sum_fit = sumlv(fitted_photons.begin(), fitted_photons.end());
        missing -= sum_raw;
        raw_n.at(comb.size() - MinNGamma()).Fill(sum_raw.M());
        fit_n.at(comb.size() - MinNGamma()).Fill(sum_fit.M());
        raw_mm.at(comb.size() - MinNGamma()).Fill(missing.M());
        raw_pE.at(comb.size() - MinNGamma()).Fill(proton->E - ParticleTypeDatabase::Proton.Mass());
        fit_pE.at(comb.size() - MinNGamma()).Fill(fitted_proton->E - ParticleTypeDatabase::Proton.Mass());
        if (fitted_proton->Candidate->VetoEnergy && fitted_proton->Candidate->Detector & Detector_t::Type_t::CB) {
            dEvE.at(comb.size() - MinNGamma()).at(fitted_proton->Candidate->FindVetoCluster()->CentralElement)
                    .Fill(fitted_proton->E - ParticleTypeDatabase::Proton.Mass(), fitted_proton->Candidate->VetoEnergy);
            dEvE_all.at(comb.size() - MinNGamma()).Fill(fitted_proton->E - ParticleTypeDatabase::Proton.Mass(),
                                                        fitted_proton->Candidate->VetoEnergy);
            dEvE_combined.at(fitted_proton->Candidate->FindVetoCluster()->CentralElement)
                    .Fill(fitted_proton->E - ParticleTypeDatabase::Proton.Mass(), fitted_proton->Candidate->VetoEnergy);
            dEvE_all_combined.Fill(fitted_proton->E - ParticleTypeDatabase::Proton.Mass(), fitted_proton->Candidate->VetoEnergy);
        }

        for (auto c = utils::makeCombination(photons, 2); !c.done(); ++c) {
            const LorentzVec sum = sumlv(c.begin(), c.end());
            raw_2.at(comb.size() - MinNGamma()).Fill(sum.M());
        }
        for (auto c = utils::makeCombination(fitted_photons, 2); !c.done(); ++c) {
            const LorentzVec sum = sumlv(c.begin(), c.end());
            fit_2.at(comb.size() - MinNGamma()).Fill(sum.M());
        }
    }
}

bool IMCombFitPlots::find_best_comb(const TTaggerHit& taggerhit, TCandidatePtrList& comb,
                                    TParticleList& fitted_photons, TParticlePtr& fitted_proton)
{
    double best_prob_fit = -std_ext::inf;
    size_t best_comb_fit = comb.size();
    TLorentzVector eta;
    TParticlePtr proton;
    TParticleList photons;

    /* kinematical checks to reduce computing time */
    const interval<double> coplanarity({-25, 25});
    const interval<double> mm = ParticleTypeDatabase::Proton.GetWindow(300);

    /* test all different combinations to find the best proton candidate */
    size_t i = 0;
    do {
        photons.clear();
        proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());  // always assume last particle is the proton
        eta.SetXYZT(0,0,0,0);
        for (size_t j = 0; j < comb.size()-1; j++) {
            photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, comb.at(j)));
            eta += TParticle(ParticleTypeDatabase::Photon, comb.at(j));
        }

        const double copl = std_ext::radian_to_degree(abs(eta.Phi() - proton->Phi())) - 180.;
        if (!coplanarity.Contains(copl))
            continue;

        LorentzVec missing = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
        missing -= eta;
        if (!mm.Contains(missing.M()))
            continue;

        /* now start with the kinematic fitting */
        auto& fit = kinfit.at(photons.size()-MinNGamma());  // choose the fitter for the right amount of photons

        auto kinfit_result = fit.DoFit(taggerhit.PhotonEnergy, proton, photons);

        if (kinfit_result.Status != APLCON::Result_Status_t::Success)
            continue;

        if (PROBABILITY_CUT)
            if (kinfit_result.Probability < PROBABILITY)
                continue;

        if (!std_ext::copy_if_greater(best_prob_fit, kinfit_result.Probability))
            continue;

        best_comb_fit = i;
        fitted_photons = fit.GetFittedPhotons();
        fitted_proton = fit.GetFittedProton();
    } while (shift_right(comb) && ++i < comb.size());

    // check if a valid combination was found
    if (best_comb_fit >= comb.size() || !isfinite(best_prob_fit))
        return false;

    // restore combinations with best probability
    while (best_comb_fit-- > 0)
        shift_right(comb);

    return true;
}

void IMCombFitPlots::Finish()
{
    int channel = 0;
    for (auto& hist : dEvE_combined) {
        TH1* h = hist.prompt->ProjectionY("_py", FIRST, LAST);
        int bins = h->GetXaxis()->GetNbins();
        for (int i = 1; i <= bins; i++)
            projections->Fill(h->GetBinCenter(i), channel, h->GetBinContent(i));
        channel++;
        delete h;
    }
}

void IMCombFitPlots::ShowResult()
{
    canvas cr2(GetName() + " raw 2comb");
    canvas cf2(GetName() + " fit 2comb");
    canvas crn(GetName() + " raw n sum");
    canvas cfn(GetName() + " fit n sum");

    for (auto& h : raw_2)
        cr2 << h.subtracted;
    for (auto& h : fit_2)
        cf2 << h.subtracted;
    for (auto& h : raw_n)
        crn << h.subtracted;
    for (auto& h : fit_n)
        cfn << h.subtracted;

    cr2 << endc;
    cf2 << endc;
    crn << endc;
    cfn << endc;
}


AUTO_REGISTER_PHYSICS(IMCombFitPlots)
