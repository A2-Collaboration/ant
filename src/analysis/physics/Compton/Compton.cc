// --------------------------- Important ---------------------------
// In Compton folder is a README file. Please read it if you
// want to use this code.

#include "Compton.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

// Defining the contructor for the Compton class
Compton::Compton(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    tagger(ExpConfig::Setup::GetDetector<TaggerDetector_t>()),
    nchannels(tagger->GetNChannels())
{
    // Bins used in histograms
    const BinSettings time_bins(2000, -200, 200);
    const BinSettings mass_bins(250, 800, 1300);
    const BinSettings angle_bins(18, 0, 360);
    const BinSettings taggerchannel_bins(nchannels);

// ------------ Histograms Created Here but not Filled ------------

//    h_WeightedTaggerTime = HistFac.makeTH1D("Weighted Tagger Time",
//                                    "t [ns]","#",
//                                    time_bins,
//                                    "h_WeightedTaggerTime"
//                                    );
//    h_MM = HistFac.makeTH1D("All particles, No Weights",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM"
//                                     );
//    h_MM1 = HistFac.makeTH1D("All Particles",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM1"
//                                     );
//    h_MM11 = HistFac.makeTH1D("Uncharged",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM11"
//                                     );
//    h_MM101 = HistFac.makeTH1D("1 Particle",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM101"
//                                     );
    h_MM111 = HistFac.makeTH1D("1 Particle, Uncharged",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MM111"
                                     );
    h3D_MM111 = HistFac.makeTH3D("1 Particle, Uncharged",
                                 "Missing Mass [MeV]",
                                 "Angle [deg]",
                                 "Tagger Channel",
                                 mass_bins,
                                 angle_bins,
                                 taggerchannel_bins ,
                                 "h3D_MM111"
                                 );
    h3D_MM111_projX = HistFac.makeTH1D("1 Particle, Uncharged",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h3D_MM111_projX"
                                     );
//    h_MM102 = HistFac.makeTH1D("2 Particles",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM102"
//                                     );
//    h_MM112 = HistFac.makeTH1D("2 Particles, "
//                                     "One is Uncharged",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM112"
//                                     );
//    h_MM1021 = HistFac.makeTH1D("2 Particles, "
//                                     "Closer Missing Mass",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM1021"
//                                     );
//    h_MM10201 = HistFac.makeTH1D("2 Particles, Coplanar",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM10201"
//                                     );
//    h_MM11201 = HistFac.makeTH1D("2 Particles, Coplanar, "
//                                     "One is Uncharged",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM11201"
//                                     );
//    h_MM10211 = HistFac.makeTH1D("2 Particles, Coplanar, "
//                                     "Closer Missing Mass",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM10211"
//                                     );
//    h_MM102001 = HistFac.makeTH1D("2 Particles, Open Ang < 15",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM102001"
//                                     );
//    h_MM112001 = HistFac.makeTH1D("2 Particles, Open Ang < 15, "
//                                     "Uncharged",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM112001"
//                                     );
//    h_MM112001_switch = HistFac.makeTH1D("2 Particles, One is Uncharged, "
//                                     "Open Ang < 15",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM112001_switch"
//                                     );
//    h_MM102011 = HistFac.makeTH1D("2 Particles, Open Ang < 15, "
//                                     "Coplanar",
//                                     "mass [MeV/c^2]","#",
//                                     mass_bins,
//                                     "h_MM102011"
//                                     );
    h_MM112011 = HistFac.makeTH1D("2 Particles, Open Ang < 15, "
                                  "Coplanar, Uncharged",
                                  "mass [MeV/c^2]","#",
                                  mass_bins,
                                  "h_MM112011"
                                     );
    h3D_MM112011 = HistFac.makeTH3D("2 Particles, Open Ang < 15, "
                                    "Coplanar, Uncharged",
                                    "Missing Mass [MeV]",
                                    "Angle [deg]",
                                    "Tagger Channel",
                                    mass_bins,
                                    angle_bins,
                                    taggerchannel_bins ,
                                    "h3D_MM112011"
                                    );
    h3D_MM112011_projX = HistFac.makeTH1D("2 Particles, Open Ang < 15, "
                                  "Coplanar, Uncharged",
                                  "mass [MeV/c^2]","#",
                                  mass_bins,
                                  "h3D_MM112011_projX"
                                     );
    h_MM112011_switch = HistFac.makeTH1D("2 Particles, One is Uncharged, "
                                         "Open Ang < 15, Coplanar",
                                         "mass [MeV/c^2]","#",
                                         mass_bins,
                                         "h_MM112011_switch"
                                         );
    h3D_MM112011_switch = HistFac.makeTH3D("2 Particles, One is Uncharged, "
                                           "Open Ang < 15, Coplanar",
                                           "Missing Mass [MeV]",
                                           "Angle [deg]",
                                           "Tagger Channel",
                                           mass_bins,
                                           angle_bins,
                                           taggerchannel_bins ,
                                           "h3D_MM112011_switch"
                                           );
    h3D_MM112011_switch_projX = HistFac.makeTH1D("2 Particles, One is Uncharged, "
                                         "Open Ang < 15, Coplanar",
                                         "mass [MeV/c^2]","#",
                                         mass_bins,
                                         "h3D_MM112011_switch_projX"
                                         );

    h_ScalarCounts = HistFac.makeTH1D("Total Counts in Tagger",
                                      "Tagger Channel","#",
                                      taggerchannel_bins,
                                      "h_ScalarCounts");

//  ---------------- Get Variables at Command Line ----------------

    // Getting the prompt random windows
    if (opts->HasOption("PR"))
        PR_windows = opts->Get<std::string>("PR","-200,7,9,19,21,200");

    // Turning sting into a vector with string entries
    const auto& PR_string_vec = std_ext::tokenize_string(PR_windows,",");

    // Specifying prompt and random ranges (see README for more info)
    promptrandom.AddRandomRange
            ({ stod(PR_string_vec.at(0)), stod(PR_string_vec.at(1)) });
    promptrandom.AddPromptRange
            ({ stod(PR_string_vec.at(2)), stod(PR_string_vec.at(3)) });
    promptrandom.AddRandomRange
            ({ stod(PR_string_vec.at(4)), stod(PR_string_vec.at(5)) });

    // Get variables at command line. The range of taggerhit energies
    // that one would like to use can be specified
    if (opts->HasOption("low"))
        tagger_energy_low = opts->Get<double>("low", 0);
    if (opts->HasOption("high"))
        tagger_energy_high = opts->Get<double>("high", 2000);

// ------------------ Getting Tagger Scalars ------------------

    slowcontrol::Variables::TaggerScalers->Request();
}

// ------------ Funtions specific to the Compton class ------------

// Checks if veto_energy (energy depositied in PID and Veto wall) meets
// threshold for particle being considered charged. Returns true if
// particle is charged and false if it is not
bool Compton::IsParticleCharged(double veto_energy)
{
    if (veto_energy < .2) { return false; }
    else { return true; }
}

// For 2 particle events.
// Checks if two particles are a photon and proton
// based on their veto energy. Returns 0 if they are
// not, returns 1 if the front is a photon and returns
// 2 is the back if a photon.
int Compton::IsChargedUncharged(const TCandidateList& candidates)
{
    // Default is that both particles are charged
    bool isfrontcharged = true;
    bool isbackcharged = true;

    // Must be 2 particles
    if (candidates.size() != 2)
    {
        // Condition not met
        LOG(ERROR) << "Size of candidates should be 2";
        return 0;
    }
    else
    {
        // Checking if particles are uncharged
        if (candidates.front().VetoEnergy < .2)
        {
            isfrontcharged = false;
        }

        if (candidates.back().VetoEnergy < .2)
        {
            isbackcharged = false;
        }

        // Return 1 or 2 if one particle is charged and the other is
        // uncharged (condition met)
        if ((isfrontcharged == false) & (isbackcharged == true))
        {
            return 1;
        }
        else if ((isfrontcharged == true) & (isbackcharged == false))
        {
            return 2;
        }

        // Condition not met
        else { return 0; }
    }
}


// Input: a candidate and the 4 momentum vectors the the
// incoming photon and proton target. Output: the missing
// mass
double Compton::GetMissingMass(const TCandidate& candidate,
                               const LorentzVec target,
                               const LorentzVec incoming)
{
    vec3 unit_vec = vec3(candidate);
    LorentzVec scattered = LorentzVec({unit_vec.x*candidate.CaloEnergy,
                                       unit_vec.y*candidate.CaloEnergy,
                                       unit_vec.z*candidate.CaloEnergy},
                                      candidate.CaloEnergy);

    // Calculating the mass of the recoil proton (missing
    // mass) from the 4 momentum vector using .M()
    // Should be 938MeV if there was a Compton event
    // involving these 2 photons
    return (incoming + target - scattered).M();
}

// Used for a 2 particle events.
// Calculates the missing mass using both particles, the outputs
// the one that is closer to the mass of a proton (938 MeV)
double Compton::GetCloserMM(const TCandidateList& candidates,
                            const LorentzVec target,
                            const LorentzVec incoming)
{
    if (candidates.size() != 2)
    {
        LOG(ERROR) << "Size of candidates should be 2";
    }
    double front_missing_mass =
            GetMissingMass(candidates.front(), target, incoming);
    double back_missing_mass =
            GetMissingMass(candidates.back(), target, incoming);


    if (abs(front_missing_mass - proton_mass) <
            abs(back_missing_mass - proton_mass))
    {
        return front_missing_mass;
    }
    else if (abs(front_missing_mass - proton_mass) >
             abs(back_missing_mass - proton_mass))
    {
        return back_missing_mass;
    }
    else
    {
        LOG(WARNING) << "Missing Masses are the same";
        return front_missing_mass;
    }
}

// For 2 particle events.
// Checks if the 2 particles are coplanar. Returns true
// if they are and false if they are not.
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

        // Default value (will not make condition)
        double diff = 0.0;

        if (front_phi > back_phi)
        {
            diff = front_phi - back_phi;
        }
        if (back_phi > front_phi)
        {
            diff = back_phi - front_phi;
        }

        // Condition: phi angles must be 180±10 deg apart. If
        // condition is met -> particles are coplanar
        if (((180*diff/M_PI) > 170) & ((180*diff/M_PI) < 190))
        {
            return true;
        }
        else { return false; }
    }
}

// For 2 particle events.
// Check the angle between the calculated and detected
// recoil proton (opening angle), should be less than
// 15 deg if event is Compton.
int Compton::IsOpeningAngle(const TCandidateList& candidates,
                            const LorentzVec target,
                            const LorentzVec incoming)
{
    // Lorentz vectors for the scattering photon (scattered)
    // and what should be the recoil proton (missing) for the
    // first and second particle in the event.
    LorentzVec front_scattered;
    LorentzVec front_missing;
    LorentzVec back_scattered;
    LorentzVec back_missing;

    vec3 front_unit_vec = vec3(candidates.front());
    vec3 back_unit_vec = vec3(candidates.back());

    front_scattered = LorentzVec({front_unit_vec.x*candidates.front().CaloEnergy,
                                  front_unit_vec.y*candidates.front().CaloEnergy,
                                  front_unit_vec.z*candidates.front().CaloEnergy},
                                  candidates.front().CaloEnergy);

    front_missing = incoming + target - front_scattered;

    back_scattered = LorentzVec({back_unit_vec.x*candidates.back().CaloEnergy,
                                 back_unit_vec.y*candidates.back().CaloEnergy,
                                 back_unit_vec.z*candidates.back().CaloEnergy},
                                 candidates.back().CaloEnergy);

    back_missing = incoming + target - back_scattered;

    // Don't know which particle is the photon, so calculate the
    // opening angle twice assuming both are the photon -> take
    // the smaller opening angle
    // THIS MIGHT NOT BE THE BEST WAY TO DO THIS
    double open_ang2 = 180*front_scattered.Angle(back_missing)/M_PI;
    double open_ang1 = 180*back_scattered.Angle(front_missing)/M_PI;

    if (open_ang1 < open_ang2)
    {
        if (open_ang1 < 15.0) { return 1; }
        else { return 0; }
    }

    else
    {
        if (open_ang2 < 15.0) { return 2; }
        else { return 0; }
    }
}

bool Compton::IsOpeningAngle2(const TCandidateList& candidates,
                     const LorentzVec target,
                     const LorentzVec incoming,
                     const int IsChargedUncharged_output)
{
    LorentzVec scattered;
    LorentzVec missing;
    LorentzVec recoil;

    vec3 front_unit_vec = vec3(candidates.front());
    vec3 back_unit_vec = vec3(candidates.back());

    double opening_angle;

    if (IsChargedUncharged_output == 1)
    {
        scattered = LorentzVec({front_unit_vec.x*candidates.front().CaloEnergy,
                                    front_unit_vec.y*candidates.front().CaloEnergy,
                                    front_unit_vec.z*candidates.front().CaloEnergy},
                                    candidates.front().CaloEnergy);

        missing = incoming + target - scattered;

        recoil = LorentzVec({back_unit_vec.x*candidates.back().CaloEnergy,
                                 back_unit_vec.y*candidates.back().CaloEnergy,
                                 back_unit_vec.z*candidates.back().CaloEnergy},
                                 candidates.back().CaloEnergy);

        opening_angle = 180*recoil.Angle(missing)/M_PI;

        if (opening_angle < 15.0 ) { return true; }

        else { return false; }
    }

    if (IsChargedUncharged_output == 2)
    {
        scattered = LorentzVec({back_unit_vec.x*candidates.back().CaloEnergy,
                                    back_unit_vec.y*candidates.back().CaloEnergy,
                                    back_unit_vec.z*candidates.back().CaloEnergy},
                                    candidates.back().CaloEnergy);

        missing = incoming + target - scattered;

        recoil = LorentzVec({front_unit_vec.x*candidates.front().CaloEnergy,
                                 front_unit_vec.y*candidates.front().CaloEnergy,
                                 front_unit_vec.z*candidates.front().CaloEnergy},
                                 candidates.front().CaloEnergy);

        opening_angle = 180*recoil.Angle(missing)/M_PI;

        if (opening_angle < 15.0 ) { return true; }

        else { return false; }
    }

    else
    {
        LOG(ERROR) << "Invalid IsChargedUncharged_output, "
                      "function returns false";
        return false;
    }
}

// ------------------------- Other Methods -------------------------

void Compton::PlotCounts()
{
    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        const auto counts = slowcontrol::Variables::TaggerScalers->
                            GetCounts().at(ch);
        h_ScalarCounts->Fill(ch,counts);
    }
}

// ----------------------- Where the Physics Happens -----------------------

void Compton::ProcessEvent(const TEvent& event, manager_t&)
{

//     --------------------- Prompt Random Stuff ---------------------

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
        // 0 and 1 in random (for calculus reasons)
        const double weight = promptrandom.FillWeight();
//        h_WeightedTaggerTime->Fill(taggerhit.Time, weight);

//     ---------- The Bulk (cuts, and filling histograms) ----------

        // Filling in the momentum 4 vec for the incoming photon
        // using tagger info
        incoming_vec = LorentzVec({0.0,0.0,taggerhit.PhotonEnergy},
                                      taggerhit.PhotonEnergy);

//            --------------- All Events ---------------

        // Looping over the candidates
        // (particles in CB) in each event
//        for (const auto& candidate : event.Reconstructed().Candidates) {

//            missing_mass = GetMissingMass(candidate, target_vec, incoming_vec);

//            // All particles, no weights
//            h_MM->Fill(missing_mass);
//            // All particles, with weights
//            h_MM1->Fill(missing_mass, weight);


//            if (Compton::IsParticleCharged(candidate.VetoEnergy) == false)
//            {
//                // All charged particles cut
//                h_MM11->Fill(missing_mass, weight);
//            }
//        }

//             -------------- 1 Particle Events --------------

        if (event.Reconstructed().Candidates.size() == 1)
        {
            for (const auto& candidate : event.Reconstructed().Candidates)
            {
                missing_mass = GetMissingMass(candidate, target_vec, incoming_vec);

                // 1 particle in event
//                h_MM101->Fill(missing_mass, weight);

                if (Compton::IsParticleCharged(candidate.VetoEnergy) == false)
                {
                    // 1 particle in event, particle is uncharged
                    h_MM111->Fill(missing_mass, weight);

                    // Filling 3D Plot
                    h3D_MM111->Fill(missing_mass, candidate.Theta,
                                    taggerhit.Channel, weight);
                }
            }
        }

//             -------------- 2 Particle Events --------------

        if (event.Reconstructed().Candidates.size() == 2)
        {
            const auto& candidates = event.Reconstructed().Candidates;

            // Using both particles to calc missing mass
//            for (const auto& candidate : candidates)
//            {
//                missing_mass = GetMissingMass(candidate, target_vec, incoming_vec);

                // 2 particles in event
//                h_MM102->Fill(missing_mass, weight);
//            }

            // Keeping only the 2 particles in which one is charged and
            // the other is not. Using the uncharged particle to calc
            // the missing mass
            if (IsChargedUncharged(candidates) == 1)
            {
//                missing_mass = GetMissingMass(candidates.front(),
//                                              target_vec, incoming_vec);

                // 2 particles in event, one is charged and the other is not
//                h_MM112->Fill(missing_mass, weight);

                // Opening angle cut
                if( IsOpeningAngle2(candidates,target_vec,incoming_vec, 1) == true )
                {
//                    missing_mass = GetMissingMass(candidates.front(),
//                                                  target_vec, incoming_vec);

                    // 2 particles in event, one is charged and the other is not,
                    // opening angle < 15 degrees
//                    h_MM112001_switch->Fill(missing_mass, weight);

                    // Colpanar cut
                    if ( IsCoplanar(candidates) == true )
                    {
                        missing_mass = GetMissingMass(candidates.front(),
                                                      target_vec, incoming_vec);

                        // 2 particles in event, one is charged and the other is not,
                        // opening angle < 15 degrees, coplanar
                        h_MM112011_switch->Fill(missing_mass, weight);

                        // Fill 3D histogram
                        h3D_MM112011_switch->Fill(missing_mass,
                                                  candidates.front().Theta,
                                                  taggerhit.Channel,
                                                  weight);
                    }
                }
            }

            // All again but with other configuration
            if (IsChargedUncharged(candidates) == 2)
            {
//                missing_mass = GetMissingMass(candidates.back(),
//                                              target_vec, incoming_vec);

                // 2 particles in event, one is charged and the other is not
//                h_MM112->Fill(missing_mass, weight);

                // Opening angle cut
                if( IsOpeningAngle2(candidates,target_vec,incoming_vec, 2) == true )
                {
//                    missing_mass = GetMissingMass(candidates.back(),
//                                                  target_vec, incoming_vec);

//                    h_MM112001_switch->Fill(missing_mass, weight);

                    // Colpanar cut
                    if ( IsCoplanar(candidates) == true )
                    {
                        missing_mass = GetMissingMass(candidates.back(),
                                                      target_vec, incoming_vec);

                        // 2 particles in event, one is charged and the other is not,
                        // opening angle < 15 degrees, coplanar
                        h_MM112011_switch->Fill(missing_mass, weight);

                        // Fill 3D histogram
                        h3D_MM112011_switch->Fill(missing_mass,
                                                  candidates.front().Theta,
                                                  taggerhit.Channel,
                                                  weight);
                    }
                }

            }

            // Getting the closer missing mass
//            closer_missing_mass = GetCloserMM
//                    (candidates, target_vec, incoming_vec);

            // 2 particles in event, only closer missing mass
//            h_MM1021->Fill(closer_missing_mass, weight);

            // Check if 2 particles in event are coplanar
            if (IsCoplanar(candidates) == true )
            {
//                for (const auto& candidate : candidates)
//                {
//                    missing_mass = GetMissingMass
//                            (candidate, target_vec, incoming_vec);

                    // 2 particles in event, event is coplanar
//                    h_MM10201->Fill(missing_mass, weight);
//                }

                // Uncharged/Charged cut
//                if (IsChargedUncharged(candidates) == 1)
//                {
//                    missing_mass = GetMissingMass(candidates.front(),
//                                                  target_vec, incoming_vec);

                    // 2 particles in event, one is charged and the
                    // other is not, event is coplanar
//                    h_MM11201->Fill(missing_mass, weight);
//                }

//                if (IsChargedUncharged(candidates) == 2)
//                {
//                    missing_mass = GetMissingMass(candidates.back(),
//                                                  target_vec, incoming_vec);

                    // 2 particles in event, one is charged and the
                    // other is not, event is coplanar
//                    h_MM11201->Fill(missing_mass,weight);
//                }

                // Closer missing mass cut
//                closer_missing_mass = GetCloserMM
//                        (candidates, target_vec, incoming_vec);

                // 2 particles in event, event is coplanar, plot
                // only the closer missing mass
//                h_MM10211->Fill(closer_missing_mass, weight);
            }

            // Check if opening angle is < 15 deg
            if (IsOpeningAngle(candidates, target_vec, incoming_vec) == 1 )
            {
//                missing_mass = GetMissingMass
//                         (candidates.front(), target_vec, incoming_vec);

                // 2 particles in event, open ang < 15
//                h_MM102001->Fill(missing_mass, weight);

                // Cut charged particles
//                if (IsParticleCharged(candidates.front().VetoEnergy) == false)
//                {
//                    missing_mass = GetMissingMass
//                             (candidates.front(), target_vec, incoming_vec);

                    // 2 particles in event, open ang < 15, uncharged
//                    h_MM112001->Fill(missing_mass, weight);
//                }

                // Keep only coplanar events
                if (IsCoplanar(candidates) == true )
                {
//                    missing_mass = GetMissingMass
//                             (candidates.front(), target_vec, incoming_vec);

                    // 2 particles in event, open ang < 15, event is coplanar
//                    h_MM102011->Fill(missing_mass, weight);

                    // Cut charged particles
                    if (IsParticleCharged(candidates.front().VetoEnergy) == false)
                    {
                        missing_mass = GetMissingMass
                                 (candidates.front(), target_vec, incoming_vec);

                        // 2 particles in event, open ang < 15, uncharged,
                        // event is coplanar
                        h_MM112011->Fill(missing_mass, weight);
                        
                        // Fill 3D histogram
                        h3D_MM112011->Fill(missing_mass,
                                           candidates.front().Theta,
                                           taggerhit.Channel,
                                           weight);
                    }
                }
            }

            // All that again but for the other configuration
            if (IsOpeningAngle(candidates, target_vec, incoming_vec) == 2 )
            {
//                missing_mass = GetMissingMass
//                         (candidates.back(), target_vec, incoming_vec);

                // 2 particles in event, open ang < 15
//                h_MM102001->Fill(missing_mass, weight);

                // Cut charged particles
//                if (IsParticleCharged(candidates.back().VetoEnergy) == false)
//                {
//                    missing_mass = GetMissingMass
//                             (candidates.back(), target_vec, incoming_vec);

                    // 2 particles in event, open ang < 15, uncharged
//                    h_MM112001->Fill(missing_mass, weight);
//                }

                // Keep only coplanar events
                if (IsCoplanar(candidates) == true )
                {
//                    missing_mass = GetMissingMass
//                             (candidates.back(), target_vec, incoming_vec);

                    // 2 particles in event, open ang < 15, event is coplanar
//                    h_MM102011->Fill(missing_mass, weight);

                    // Cut charged particles
                    if (IsParticleCharged(candidates.front().VetoEnergy) == false)
                    {
                        missing_mass = GetMissingMass
                                 (candidates.back(), target_vec, incoming_vec);

                        // 2 particles in event, open ang < 15, uncharged,
                        // event is coplanar
                        h_MM112011->Fill(missing_mass, weight);

                        // Fill 3D histogram
                        h3D_MM112011->Fill(missing_mass,
                                           candidates.front().Theta,
                                           taggerhit.Channel,
                                           weight);
                    }
                }
            }
        }
    }

    if(slowcontrol::Variables::TaggerScalers->HasChanged())
    {
        PlotCounts();
    }

    h3D_MM111_projX =
            h3D_MM111->ProjectionX();
    h3D_MM112011_projX =
            h3D_MM112011->ProjectionX();
    h3D_MM112011_switch_projX =
            h3D_MM112011_switch->ProjectionX();
}

// ---------------------- Outputing the Histograms ----------------------

void Compton::ShowResult()
{

//    ant::canvas(GetName()+": Tagger Time Plots")
//            << h_WeightedTaggerTime
//            << endc; // actually draws the canvas

//    ant::canvas(GetName()+": Preliminary Cuts")
//            << h_MM
//            << h_MM1
//            << h_MM11
//            << endc;

    ant::canvas(GetName()+": 1 Particle Events")
//            << h_MM101
            << h_MM111
            << h3D_MM111
            << h3D_MM111_projX
            << endc;

//    ant::canvas(GetName()+": 2 Particle Events")
//            << h_MM102
//            << h_MM112
//            << h_MM1021
//            << endc;

//    ant::canvas(GetName()+": Coplanar Events")
//            << h_MM10201
//            << h_MM11201
//            << h_MM10211
//            << endc;

    ant::canvas(GetName()+": Events with Opening Angle < 15")
//            << h_MM102001
//            << h_MM112001
//            << h_MM102011
            << h_MM112011
            << h3D_MM112011
            << h3D_MM112011_projX
//            << h_MM112001_switch
            << h_MM112011_switch
            << h3D_MM112011_switch
            << h3D_MM112011_switch_projX
            << endc;

    ant::canvas(GetName()+": Scalar Counts")
            << h_ScalarCounts
            << endc;
}

// ---------------- Registering the Class ----------------

// A macro that registers the Compton class with Ant so that
// you can call this class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
