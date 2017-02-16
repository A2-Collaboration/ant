#include "PID_Energy.h"

#include "expconfig/ExpConfig.h"

#include "base/std_ext/vector.h"

#include "base/Logger.h"

#include <numeric>  // std::accumulate

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

template<typename T>
void PID_Energy::shift_right(std::vector<T>& v)
{
    std::rotate(v.begin(), v.end() -1, v.end());
}

APLCON::Fit_Settings_t PID_Energy::MakeFitSettings(unsigned max_iterations)
{
    auto settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = max_iterations;
    return settings;
}

PID_Energy::PID_Energy(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    useMIP(opts->Get<bool>("UseMIP", false)),
    model(make_shared<utils::UncertaintyModels::FitterSergey>()),
    kinfit("kinfit", 3, model, true, MakeFitSettings(20))
{
    promptrandom.AddPromptRange({-3, 2});

    const auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);
    const auto nChannels = detector->GetNChannels();

    const BinSettings pid_channels(nChannels);
    const BinSettings pid_rawvalues(300);
    const BinSettings energybins(500, 0, 10);


    h_pedestals = HistFac.makeTH2D(
                      "PID Pedestals",
                      "Raw ADC value",
                      "#",
                      pid_rawvalues,
                      pid_channels,
                      "Pedestals");

    h_bananas =
            HistFac.makeTH3D(
                "PID Bananas",
                "CB Energy / MeV",
                "PID Energy / MeV",
                "Channel",
                BinSettings(400,0,800),
                BinSettings(200,0,18),
                pid_channels,
                "Bananas"
                );

    h_mip = HistFac.makeTH2D(
                "PID Minimum Ionizing Peak",
                "PID Energy / MeV",
                "Channel",
                energybins,
                pid_channels,
                "MIP");

    for(unsigned ch=0;ch<nChannels;ch++) {
        stringstream ss;
        ss << "Ch" << ch;
        h_perChannel.push_back(
                    PerChannel_t(HistogramFactory(ss.str(), HistFac, ss.str())));
    }

    kinfit.SetZVertexSigma(3);

    if (useMIP)
        LOG(INFO) << "Create PID Calibration histograms for Minimum Ionizing Peak method";
}

PID_Energy::PerChannel_t::PerChannel_t(HistogramFactory HistFac)
{
    const BinSettings cb_energy(400,0,800);
    const BinSettings pid_timing(300,-300,700);
    const BinSettings pid_rawvalues(300);
    const BinSettings pid_energy(150,0,30);

    PedestalTiming = HistFac.makeTH2D(
                         "PID Pedestals Timing",
                         "Timing / ns",
                         "Raw ADC value",
                         pid_timing,
                         pid_rawvalues,
                         "PedestalTiming");

    PedestalNoTiming = HistFac.makeTH1D(
                           "PID Pedestals No Timing",
                           "Raw ADC value",
                           "#",
                           pid_rawvalues,
                           "PedestalNoTiming");

    Banana = HistFac.makeTH2D(
                 "PID Banana",
                 "CB Energy / MeV",
                 "PID Energy / MeV",
                 cb_energy,
                 pid_energy,
                 "Banana"
                 );

    BananaRaw = HistFac.makeTH2D(
                    "PID Banana Raw",
                    "CB Energy / MeV",
                    "PID ADC Value",
                    cb_energy,
                    BinSettings(300,0,2000),
                    "BananaRaw"
                    );


    TDCMultiplicity = HistFac.makeTH1D("PID TDC Multiplicity", "nHits", "#", BinSettings(10), "TDCMultiplicity");
    QDCMultiplicity = HistFac.makeTH1D("PID QDC Multiplicity", "nHits", "#", BinSettings(10), "QDCMultiplicity");
}

void PID_Energy::ProcessEvent(const TEvent& event, manager_t&)
{
    // pedestals, best determined from clusters with energy information only

    struct hitmapping_t {
        vector<double> Pedestals;
        vector<double> Timings;
    };

    std::map<unsigned, hitmapping_t> hits;

    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        if(readhit.DetectorType != Detector_t::Type_t::PID)
            continue;

        auto& item = hits[readhit.Channel];

        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            std_ext::concatenate(item.Pedestals, readhit.Converted);
        }
        else if(readhit.ChannelType == Channel_t::Type_t::Timing) {
            std_ext::concatenate(item.Timings, readhit.Values); // passed the timing window!
        }
    }

    for(const auto& it_hit : hits) {

        const auto channel = it_hit.first;
        const hitmapping_t& item = it_hit.second;

        PerChannel_t& h = h_perChannel[channel];

        h.QDCMultiplicity->Fill(item.Pedestals.size());
        h.TDCMultiplicity->Fill(item.Timings.size());


        if(item.Pedestals.size() != 1)
            continue;
        if(item.Timings.size()>1)
            continue;

        const auto& pedestal = item.Pedestals.front();

        h_pedestals->Fill(pedestal, channel);

        if(item.Timings.size()==1)
            h.PedestalTiming->Fill(item.Timings.front(), pedestal);
        else
            h.PedestalNoTiming->Fill(pedestal);

    }

    // bananas per channel histograms
    for(const auto& candidate : event.Reconstructed().Candidates) {
        // only candidates with one cluster in CB and one cluster in PID
        if(candidate.Clusters.size() != 2)
            continue;
        const bool cb_and_pid = candidate.Detector & Detector_t::Type_t::CB &&
                                candidate.Detector & Detector_t::Type_t::PID;
        if(!cb_and_pid)
            continue;

        // search for PID cluster
        const auto& pid_cluster = candidate.FindFirstCluster(Detector_t::Type_t::PID);

        h_bananas->Fill(candidate.CaloEnergy,
                        candidate.VetoEnergy,
                        pid_cluster->CentralElement);

        // per channel histograms
        PerChannel_t& h = h_perChannel[pid_cluster->CentralElement];

        // fill the banana
        h.Banana->Fill(candidate.CaloEnergy,
                       candidate.VetoEnergy);

        // is there an pedestal available?
        const auto it_hit = hits.find(pid_cluster->CentralElement);
        if(it_hit == hits.end()) {
            continue;
        }
        const auto& pedestals = it_hit->second.Pedestals;
        if(pedestals.size() != 1)
            continue;

        const auto& pedestal = pedestals.front();

        h.BananaRaw->Fill(candidate.CaloEnergy, pedestal);
        //h.BananaTiming->Fill(candidate.ClusterEnergy(), candidate.VetoEnergy(), timing);

    }


    if (useMIP)
        ProcessMIP(event);
}

void PID_Energy::ProcessMIP(const TEvent& event)
{
    // analyze e+ e- gamma events, determine proton via kinematic fit
    const auto& data = event.Reconstructed();
    const auto& cands = data.Candidates;
    const bool MC = data.ID.isSet(TID::Flags_t::MC);

    if (cands.size() != 4)
        return;

    if (MC)
        if (data.Trigger.CBEnergySum <= 550)
            return;

    TParticlePtr proton;
    TCandidatePtrList comb;
    for (auto p : cands.get_iter())
        comb.emplace_back(p);

    // require at least 2 candidates with PID/Veto entries
    if (std::count_if(comb.begin(), comb.end(), [](TCandidatePtr c){ return c->VetoEnergy; }) < 2)
        return;

    double CBAvgTime = event.Reconstructed().Trigger.CBTiming;
    if (MC)
        CBAvgTime = 0;
    if (!isfinite(CBAvgTime))
        return;

    TParticleList photons;
    TParticleList fitted_photons;  // used to store the kinfitted photon information
    double best_prob_fit = -std_ext::inf;
    size_t best_comb_fit = cands.size();
    for (const TTaggerHit& taggerhit : data.TaggerHits) {  // loop over all tagger hits
        promptrandom.SetTaggerHit(taggerhit.Time - CBAvgTime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        for (size_t i = 0; i < cands.size(); i++) {  // loop to test all different combinations
            // ensure the possible proton candidate is kinematically allowed
            if (std_ext::radian_to_degree(comb.back()->Theta) > 90.) {
                shift_right(comb);
                continue;
            }

            // require 2 PID entries for the meson candidate
            if (std::count_if(comb.begin(), comb.end()-1, [](TCandidatePtr c){ return c->VetoEnergy; }) < 2) {
                shift_right(comb);
                continue;
            }

            photons.clear();
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, comb.back());
            for (size_t j = 0; j < comb.size()-1; j++)
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, comb.at(j)));

            // do the fitting and check if the combination is better than the previous best
            if (!doFit_checkProb(taggerhit, proton, photons, best_prob_fit, fitted_photons)) {
                shift_right(comb);
                continue;
            }

            best_comb_fit = i;

            shift_right(comb);
        }
    }

    if (best_comb_fit >= cands.size() || !isfinite(best_prob_fit))
        return;

    // cut on the fit probability
    if (best_prob_fit < .05)
        return;

    // restore combinations with best probability
    while (best_comb_fit-- > 0)
        shift_right(comb);

    // sort the meson final state according to their Veto energies
    sort(comb.begin(), comb.end()-1,
         [] (const TCandidatePtr& a, const TCandidatePtr& b) {
            return a->VetoEnergy > b->VetoEnergy;
         });

    const TCandidatePtr& l1 = comb.at(0);
    const TCandidatePtr& l2 = comb.at(1);
    // suppress conversion decays
    if (l1->FindVetoCluster()->CentralElement == l2->FindVetoCluster()->CentralElement)
        return;

    // calculate IM of fitted photons for best combination
    TLorentzVector im;
    im = accumulate(fitted_photons.begin(), fitted_photons.end(), TLorentzVector(0,0,0,0),
                    [](TLorentzVector sum, TParticlePtr p){ return sum += *p; });

    // cut on IM(e+e-g) to suppress charged pions and other unwanted higher energetic decay channels
    if (im.M() > 600)
        return;

    // sort fitted photons according to their Veto energies
    sort(fitted_photons.begin(), fitted_photons.end()-1,
         [] (const TParticlePtr& a, const TParticlePtr& b) {
            return a->Candidate->VetoEnergy > b->Candidate->VetoEnergy;
         });

    // cut on IM(e+e-) to suppress charged pions
    im = *(fitted_photons.at(0)) + *(fitted_photons.at(1));
    if (im.M() > 300)
        return;

    // fill calibration histogram
    for (const TCandidatePtr& c : comb)
        if (c->VetoEnergy && c->Detector & Detector_t::Type_t::CB)
            h_mip->Fill(c->VetoEnergy, c->FindVetoCluster()->CentralElement);
}

void PID_Energy::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << h_pedestals
            << endc;
    canvas c_bananas(GetName()+": Bananas");
    c_bananas << drawoption("colz");
    for(auto& h : h_perChannel)
        c_bananas << h.Banana;
    c_bananas << endc;

    if (useMIP)
        canvas(GetName()+": MIP") << drawoption("colz") << h_mip << endc;
}

bool PID_Energy::doFit_checkProb(const TTaggerHit& taggerhit,
                                 const TParticlePtr proton,
                                 const TParticleList photons,
                                 double& best_prob_fit,
                                 TParticleList& fit_photons)
{
    TLorentzVector meson(0,0,0,0);

    for (const auto& g : photons)
        meson += *g;

    /* kinematical checks to reduce computing time */
    const interval<double> coplanarity({-25, 25});
    const interval<double> mm = ParticleTypeDatabase::Proton.GetWindow(300);

    const double copl = std_ext::radian_to_degree(abs(meson.Phi() - proton->Phi())) - 180.;
    if (!coplanarity.Contains(copl))
        return false;

    LorentzVec missing = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
    missing -= meson;
    if (!mm.Contains(missing.M()))
        return false;


    // now start with the kinematic fitting
    kinfit.SetEgammaBeam(taggerhit.PhotonEnergy);
    kinfit.SetProton(proton);
    kinfit.SetPhotons(photons);

    auto kinfit_result = kinfit.DoFit();

    if (kinfit_result.Status != APLCON::Result_Status_t::Success)
        return false;

    const double prob = kinfit_result.Probability;


    if (!std_ext::copy_if_greater(best_prob_fit, prob))
        return false;

    // get the fitted photon information
    fit_photons.clear();
    fit_photons = kinfit.GetFittedPhotons();

    return true;
}

AUTO_REGISTER_PHYSICS(PID_Energy)
