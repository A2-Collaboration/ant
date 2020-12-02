/*
 * Code to go in final Compton class (but I don't want to run every time)
 * /////////////////////////////////////////////////////////////////////
*/

// In consructor, creating histogram
h_TaggerCBSubtaction = HistFac.makeTH1D("Tagger CB Subtaction",
                                "t [ns]","#",
                                time_bins,
                                "h_TaggerCBSubtaction"
                                );

// In ProcessEvent, plot you would use to get Prompt and
// random time windows (to figure out later)
h_TaggerCBSubtaction->Fill(CorrectedTaggerTime);

/*
 * Code I might use or is helpful for other reasons
 * ////////////////////////////////////////////////
*/

// In h file
double GetMissingMass(const double& incoming_ph_energy,
                              const double& scattered_ph_energy,
                              const double& theta);

// In constructor, creating a histogram
h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
                                "t [ns]","#",     // xlabel, ylabel
                                time_bins, // our binnings
                                "h_TaggerTime"    // ROOT object name, auto-generated if omitted
                                );

// Trying to get missing mass using the formula
// Didn't work, but maybe useful to revisit
double Compton::GetMissingMass(const double& incoming_ph_energy,
                            const double& scattered_ph_energy,
                            const double& theta) {
    // Calculates missing mass given the incoming photon energy,
    // scattered photon energy, and scattered photon angle
    return (incoming_ph_energy * scattered_ph_energy)*(1.0 - cos(theta))
            /(incoming_ph_energy - scattered_ph_energy);
}

// In Candidate for loop. Calculating missing mass
double MissingMass = GetMissingMass(taggerhit.PhotonEnergy,
                               candidate.CaloEnergy, candidate.Theta);
h_MissingMass->Fill(MissingMass);

if (event.Reconstructed().Candidates.size() == 2)
{
    const auto& candidates = event.Reconstructed().Candidates;

    vec3 front_unit_vec = vec3(candidates.front());
    vec3 back_unit_vec = vec3(candidates.back());

    LorentzVec front_vec = LorentzVec({front_unit_vec.x*candidates.front().CaloEnergy,
                                  front_unit_vec.y*candidates.front().CaloEnergy,
                                  front_unit_vec.z*candidates.front().CaloEnergy},
                                 candidates.front().CaloEnergy);

    LorentzVec back_vec = LorentzVec({back_unit_vec.x*candidates.back().CaloEnergy,
                                 back_unit_vec.y*candidates.back().CaloEnergy,
                                 back_unit_vec.z*candidates.back().CaloEnergy},
                                candidates.back().CaloEnergy);

    h_Pi0->Fill((front_vec + back_vec).M());

}

// Opening Angle Filter
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
front_missing = incoming_vec + target_vec - front_scattered;

back_scattered = LorentzVec({back_unit_vec.x*candidates.back().CaloEnergy,
                             back_unit_vec.y*candidates.back().CaloEnergy,
                             back_unit_vec.z*candidates.back().CaloEnergy},
                            candidates.back().CaloEnergy);
back_missing = incoming_vec + target_vec - back_scattered;

open_ang = front_scattered.Angle(back_missing);
h_OpeningAngle->Fill(180*open_ang/M_PI);
open_ang = back_scattered.Angle(front_missing);
h_OpeningAngle->Fill(180*open_ang/M_PI);

const BinSettings mass_bins2(500, -20, 2000);
const BinSettings angle_bins(100 , 0 , 200);
const BinSettings energy_bins(200 , 0 , 1000);

double Compton::GetMissingMass(const TCandidate& candidate,
                 LorentzVec target, LorentzVec incoming)
{
    LorentzVec scattered = LorentzVec(vec3(candidate),candidate.CaloEnergy);

    const TParticle target_particle(ParticleTypeDatabase::Proton,target);
    const TParticle incoming_particle(ParticleTypeDatabase::Photon,incoming);
    const TParticle scattered_particle(ParticleTypeDatabase::Photon,scattered);

    // Calculating the mass of the recoil proton from
    // the 4 momentum vector using .M()
    // Should be 938MeV if there was a Compton
    // event involving these 2 photons
    return (incoming_particle + target_particle - scattered_particle).M();
}

//LorentzVec scattered = LorentzVec(vec3(candidate),candidate.CaloEnergy);
vec3 unit_vec = vec3(candidate);
LorentzVec scattered = LorentzVec({unit_vec.x*candidate.CaloEnergy,unit_vec.y*candidate.CaloEnergy,unit_vec.z*candidate.CaloEnergy},candidate.CaloEnergy);
h_ScatteredMass->Fill(scattered.M());

const TParticle scattered_particle(ParticleTypeDatabase::Photon,scattered);
h_ScatteredMass2->Fill(scattered_particle.M());

h_Theta->Fill(180*candidate.Theta/M_PI);


