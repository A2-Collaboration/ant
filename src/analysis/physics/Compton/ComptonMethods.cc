#include "Compton.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

// Defining the contructor for the Compton class
Compton::Compton(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    // Histograms are created here but not filled

    const BinSettings time_bins(2000, -200, 200);
    BinSettings mass_bins(1000, 700, 1100);
    const BinSettings angle_bins1(60 , 120 , 240);
    const BinSettings angle_bins2(90 , 0 , 180);

    h_PromptRandomWithTriggerSimulation = HistFac.makeTH1D("PromptRandom with TriggerSimulation",
                                    "t [ns]","#",
                                    time_bins,
                                    "h_PromptRandomWithTriggerSimulation"
                                    );
    h_MissingMass = HistFac.makeTH1D("All Taggerhits, all Candidates",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass"
                                     );
    h_MissingMass1 = HistFac.makeTH1D("Weighted Taggerhits, all Candidates",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass1"
                                     );
    h_MissingMass01 = HistFac.makeTH1D("All Taggerhits, no Veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass01"
                                     );
    h_MissingMass11 = HistFac.makeTH1D("Weighted Taggerhits, no Veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass11"
                                     );
    h_MissingMass001 = HistFac.makeTH1D("All Taggerhits, 1 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass001"
                                     );
    h_MissingMass101 = HistFac.makeTH1D("Weighted Taggerhits, 1 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass101"
                                     );
    h_MissingMass011 = HistFac.makeTH1D("All Taggerhits, 1 particle event with no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass011"
                                     );
    h_MissingMass111 = HistFac.makeTH1D("Weighted Taggerhits, 1 particle event with no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass111"
                                     );
    h_MissingMass002 = HistFac.makeTH1D("All Taggerhis, 2 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass002"
                                     );
    h_MissingMass102 = HistFac.makeTH1D("Weighted Taggerhis, 2 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass102"
                                     );
    h_MissingMass012 = HistFac.makeTH1D("All Taggerhits, 2 particle event, only one has no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass012"
                                     );
    h_MissingMass112 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle event, only one has no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass112"
                                     );
    h_MissingMass0021 = HistFac.makeTH1D("All Taggerhits, 2 particle event, only plot particle with "
                                         "closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass0021"
                                     );
    h_MissingMass1021 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle event, only plot particle with "
                                         "closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass1021"
                                     );
    h_MissingMass00201 = HistFac.makeTH1D("All Taggerhits, 2 particle coplanar event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass00201"
                                     );
    h_MissingMass10201 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle coplanar event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass10201"
                                     );
    h_MissingMass01201 = HistFac.makeTH1D("All Taggerhits, 2 particle coplanar event, no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass01201"
                                     );
    h_MissingMass11201 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle coplanar event, no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass11201"
                                     );
    h_MissingMass00211 = HistFac.makeTH1D("All Taggerhits, 2 particle coplanar event,"
                                          " plot only particle with closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass00211"
                                     );
    h_MissingMass10211 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle coplanar event,"
                                          " plot only particle with closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass10211"
                                     );

    h_CoplanarAngle = HistFac.makeTH1D("Coplanar Angle",
                                     "anlge [degrees]","#",
                                     angle_bins1,
                                     "h_CoplanarAngle"
                                     );
    h_OpeningAngle = HistFac.makeTH1D("Opening Angle",
                                     "anlge [degrees]","#",
                                     angle_bins2,
                                     "h_OpeningAngle"
                                     );


    // Prompt and random windows. Must be selected based on
    // tagger time plots **(need to figure out how to code this)**
    promptrandom.AddPromptRange({ 9, 19 }); // in nanoseconds
    promptrandom.AddRandomRange({ -200, 7 });
    promptrandom.AddRandomRange({ 21, 200 });

    // Variables that can be modified at
    // the command line
    tagger_energy_low = 0;
    tagger_energy_high = 1000;
    if (opts->HasOption("low"))
        tagger_energy_low = opts->Get<double>("low", 0);
    if (opts->HasOption("high"))
        tagger_energy_high = opts->Get<double>("high", 1000);

}
// Checks if veto_energy meets threshold for particle
// being considered charged. Returns true if particle
// is charged and false if it is not
bool Compton::IsParticleCharged(double veto_energy)
{
    if (veto_energy < .2) { return false; }
    else { return true; }
}

// Checks if two particles are a photon and proton
// based on their veto energy. Returns 0 if they are
// not, returns 1 if the front is a photon and returns
// 2 is the back if a photon.
int Compton::IsPhotonProton(const TCandidateList& candidates)
{
    bool isfrontcharged = true;
    bool isbackcharged = true;
    // Will stay a null pointer if TCandidatePtrList
    // inputed is not a photon and proton. Otherwise,
    // will be assigned the photon candidate

    // Must be 2 particles
    if (candidates.size() != 2)
    {
        // Condition not met, return 0
        LOG(WARNING) << "Size of candidates should be 2";
        return 0;
    }
    else
    {
        if (candidates.front().VetoEnergy < .2)
        {
            isfrontcharged = false;
        }

        if (candidates.back().VetoEnergy < .2)
        {
            isbackcharged = false;
        }

        if ((isfrontcharged == false) & (isbackcharged == true))
        {
            return 1;
        }
        else if ((isfrontcharged == true) & (isbackcharged == false))
        {
            return 2;
        }

        // Condition not met, photon is a nullptr
        else { return 0; }
    }
}


// Input: a candidate and the 4 momentum vectors the the
// incoming photon and proton target. Output: the missing
// mass
double Compton::GetMissingMass(const TCandidate& candidate,
                 LorentzVec target, LorentzVec incoming_ph)
{
    // Momentum 4 vector for scattered photon
    scattered_ph_vec = LorentzVec(vec3(candidate),
                                  candidate.CaloEnergy);
    // Calculating the momentum 4 vector for the possible
    // recoil proton
    recoil_pr_vec = incoming_ph + target - scattered_ph_vec;
    // Calculating the mass of the recoil proton from
    // the 4 momentum vector using .M()
    // Should be 938MeV if there was a Compton
    // event involving these 2 photons
    return recoil_pr_vec.M();
}

double Compton::GetCloserMM
(const TCandidateList& candidates, const LorentzVec target, LorentzVec incoming_ph)
{
    if (candidates.size() != 2)
    {
        LOG(ERROR) << "Size of candidates should be 2";
    }
    double front_missing_mass =
            GetMissingMass(candidates.front(), target, incoming_ph);
    double back_missing_mass =
            GetMissingMass(candidates.back(), target, incoming_ph);


    if (abs(front_missing_mass - proton_mass) < abs(back_missing_mass - proton_mass))
    {
        return front_missing_mass;
    }
    else if (abs(front_missing_mass - proton_mass) > abs(back_missing_mass - proton_mass))
    {
        return back_missing_mass;
    }
    else
    {
        LOG(WARNING) << "Missing Masses are the same";
        return front_missing_mass;
    }
}

// Coplanar filter
bool Compton::IsCoplanar(const TCandidateList& candidates)
{
    if (candidates.size() != 2)
    {
        LOG(WARNING) << "Size of candidates should be 2";
        return false;
    }

    else
    {
        double front_phi = candidates.front().Phi;
        double back_phi = candidates.back().Phi;

        double diff;

        if (front_phi > back_phi)
        {
            diff = front_phi - back_phi;
        }
        if (back_phi > front_phi)
        {
            diff = back_phi - front_phi;
        }

        // plus or minus 10 degrees from 180
        if (((180*diff/M_PI) > 170) & ((180*diff/M_PI) < 190))
        {
            return true;
        }
        else { return false; }
    }
}

void Compton::ProcessEvent(const TEvent& event, manager_t&)
{
    // Runs ProcessEvent function in TriggerSimulation file which
    // does the calculations
    triggersimu.ProcessEvent(event);

    for (const auto& taggerhit : event.Reconstructed().TaggerHits)
    {

        // Skipping taggerhits outside the specified energy range
        if ((taggerhit.PhotonEnergy < tagger_energy_low) ||
                (taggerhit.PhotonEnergy > tagger_energy_high))
        {
            continue;
        }

        // Apply trigger simulation to tagger hits
        // This subtracts a weighted time from the CB (see wiki)
        const auto& CorrectedTaggerTime =
                triggersimu.GetCorrectedTaggerTime(taggerhit);

        // This assigns weights to the TaggerHits based on which
        // time window they fall into
        promptrandom.SetTaggerTime(CorrectedTaggerTime);

        // When the Tagger hit is in neither the prompt or random
        // window, then skip
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // Plot taggerhits with weights. Weight of 1.0 in
        // prompt, weight of 0.0 in outside and weight of between
        // 0 and 1 in random
        const double weight = promptrandom.FillWeight();
        h_PromptRandomWithTriggerSimulation->Fill(taggerhit.Time, weight);

        // Calculating the momentum 4 vec for the incoming photon
        incoming_ph_vec = LorentzVec ({0,0,taggerhit.PhotonEnergy},
                                      taggerhit.PhotonEnergy);

        // Events with all numbers of particles. Looping over the
        // candidates in an event
        for (const auto& candidate : event.Reconstructed().Candidates) {

            missing_mass = GetMissingMass(candidate, target_vec, incoming_ph_vec);

            h_MissingMass->Fill(missing_mass);
            h_MissingMass1->Fill(missing_mass,weight);

            // Filter 2: Veto
            if (Compton::IsParticleCharged(candidate.VetoEnergy) == false)
            {
                h_MissingMass01->Fill(missing_mass);
                h_MissingMass11->Fill(missing_mass,weight);
            }

            // Opening angle
            //open_ang = recoil_pr_vec.Angle(scattered_ph_vec);
            //h_OpeningAngle->Fill(180*open_ang/M_PI);


        }
        // Filter 3: check if there were two particles
        // present in the event
        if (event.Reconstructed().Candidates.size() == 2)
        {
            const auto& candidates = event.Reconstructed().Candidates;
            //const auto& candidates_ptr = *candidates;

            // Using both candidates to calc missing mass
            for (const auto& candidate : candidates)
            {
                missing_mass = GetMissingMass(candidate, target_vec, incoming_ph_vec);

                // 2 particles in event, with and without weights
                h_MissingMass002->Fill(missing_mass);
                h_MissingMass102->Fill(missing_mass,weight);
            }

            // Plotting only if one is a photon and one is a proton. Only
            // the photon is used to calculate the missing mass
            if (IsPhotonProton(candidates) == 1)
            {
                missing_mass = GetMissingMass(candidates.front(),
                                              target_vec, incoming_ph_vec);
                h_MissingMass012->Fill(missing_mass);
                h_MissingMass112->Fill(missing_mass,weight);
            }

            if (IsPhotonProton(candidates) == 2)
            {
                missing_mass = GetMissingMass(candidates.back(),
                                              target_vec, incoming_ph_vec);
                h_MissingMass012->Fill(missing_mass);
                h_MissingMass112->Fill(missing_mass,weight);
            }

            // Only plotting the candidate that gives the closer
            // missing mass
            closer_missing_mass = GetCloserMM
                    (candidates, target_vec, incoming_ph_vec);
            h_MissingMass0021->Fill(closer_missing_mass);
            h_MissingMass1021->Fill(closer_missing_mass,weight);

            // Check if 2 particles in event are coplanar
            if (IsCoplanar(candidates) == true )
            {
                for (const auto& candidate : candidates)
                {
                    missing_mass = GetMissingMass
                            (candidate, target_vec, incoming_ph_vec);
                    h_MissingMass00201->Fill(missing_mass);
                    h_MissingMass10201->Fill(missing_mass,weight);
                }

                // Veto filter
                if (IsPhotonProton(candidates) == 1)
                {
                    missing_mass = GetMissingMass(candidates.front(),
                                                  target_vec, incoming_ph_vec);
                    h_MissingMass012->Fill(missing_mass);
                    h_MissingMass112->Fill(missing_mass,weight);
                }

                if (IsPhotonProton(candidates) == 2)
                {
                    missing_mass = GetMissingMass(candidates.back(),
                                                  target_vec, incoming_ph_vec);
                    h_MissingMass012->Fill(missing_mass);
                    h_MissingMass112->Fill(missing_mass,weight);
                }

                // Closer missing mass filter
                closer_missing_mass =
                        GetCloserMM
                        (candidates, target_vec, incoming_ph_vec);
                h_MissingMass00211->Fill(closer_missing_mass);
                h_MissingMass10211->Fill(closer_missing_mass,weight);
            }
        }

        if (event.Reconstructed().Candidates.size() == 1)
        {
            for (const auto& candidate : event.Reconstructed().Candidates)
            {
                missing_mass = GetMissingMass(candidate, target_vec, incoming_ph_vec);

                // 1 particle, no Veto, with and without weights
                h_MissingMass001->Fill(missing_mass);
                h_MissingMass101->Fill(missing_mass,weight);

                if (Compton::IsParticleCharged(candidate.VetoEnergy) == false)
                {
                    // 1 uncharged particle with and without weights
                    h_MissingMass011->Fill(missing_mass);
                    h_MissingMass111->Fill(missing_mass,weight);
                }
            }
        }
    }
}


void Compton::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time Plots, incoming photon energy range: "
                          )
            << h_PromptRandomWithTriggerSimulation
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Missing Mass Plots, incoming photon energy range: "
                          )
            << h_MissingMass
            << h_MissingMass1
            << h_MissingMass01
            << h_MissingMass11
            << endc;

    ant::canvas(GetName()+": Missing Mass Plots 2, incoming photon energy range: "
                          )
            << h_MissingMass001
            << h_MissingMass101
            << h_MissingMass011
            << h_MissingMass111
            << h_MissingMass002
            << h_MissingMass102
            << h_MissingMass012
            << h_MissingMass112
            << h_MissingMass0021
            << h_MissingMass1021
            << endc;

    ant::canvas(GetName()+": Coplanar Missing Mass Plots, incoming photon energy range: "
                          )
            << h_MissingMass00201
            << h_MissingMass10201
            << h_MissingMass01201
            << h_MissingMass11201
            << h_MissingMass00211
            << h_MissingMass10211
            << endc;

    ant::canvas(GetName()+": Coplaner and Opening Angle, incoming photon energy range: "
                          )
            << h_CoplanarAngle
            << h_OpeningAngle
            << endc;
}

// A macro that registers the Compton class with Ant so that
// you can call this class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
