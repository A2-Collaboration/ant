#pragma once

#include "expconfig/Detector_t.h"

#include "tree/TCluster.h"

#include "TVector3.h"
#include "TMath.h"

#include <vector>
#include <list>
#include <set>

namespace ant {
namespace reconstruct {
namespace clustering {

struct crystal_t  {
    double Energy;
    const ClusterDetector_t::Element_t* Element;
    const TClusterHit* Hit;
    crystal_t(
            double energy,
            const ClusterDetector_t::Element_t* element,
            const TClusterHit* hit
            )
        :
          Energy(energy),
          Element(element),
          Hit(hit)
    {}
};

inline bool operator< (const crystal_t& lhs, const crystal_t& rhs){
    return lhs.Energy>rhs.Energy;
}

struct bump_t {
    TVector3 Position;
    std::vector<double> Weights;
    size_t MaxIndex; // index of highest weight
};

static double calc_total_energy(const std::vector< crystal_t >& cluster) {
    double energy = 0;
    for(size_t i=0;i<cluster.size();i++) {
        energy += cluster[i].Energy;
    }
    return energy;
}

//static double opening_angle(const crystal_t& c1, const crystal_t& c2)  {
//  // use TMath::ACos since it catches some NaN cases due to double rounding issues
//  return TMath::RadToDeg()*TMath::ACos(c1.Position.Unit()*c2.Position.Unit());
//}

static double calc_energy_weight(const double energy, const double total_energy) {
    double wgtE = 4.0 + TMath::Log(energy / total_energy); /// \todo use optimal cutoff value
    return wgtE<0 ? 0 : wgtE;
}

static void calc_bump_weights(const std::vector<crystal_t>& cluster, bump_t& bump) {
    double w_sum = 0;
    for(size_t i=0;i<cluster.size();i++) {
        double r = (bump.Position - cluster[i].Element->Position).Mag();
        double w = cluster[i].Energy*TMath::Exp(-2.5*r/cluster[i].Element->MoliereRadius);
        bump.Weights[i] = w;
        w_sum += w;
    }
    // normalize weights and find index of highest weight
    // (important for merging later)
    double w_max = 0;
    size_t i_max = 0;
    for(size_t i=0;i<cluster.size();i++) {
        bump.Weights[i] /= w_sum;
        if(w_max<bump.Weights[i]) {
            i_max = i;
            w_max = bump.Weights[i];
        }
    }
    bump.MaxIndex = i_max;
}

static void update_bump_position(const std::vector<crystal_t>& cluster, bump_t& bump) {
    double bump_energy = 0;
    for(size_t i=0;i<cluster.size();i++) {
        bump_energy += bump.Weights[i] * cluster[i].Energy;
    }
    TVector3 position(0,0,0);
    double w_sum = 0;
    for(size_t i=0;i<cluster.size();i++) {
        double energy = bump.Weights[i] * cluster[i].Energy;
        double w = calc_energy_weight(energy, bump_energy);
        position += w * cluster[i].Element->Position;
        w_sum += w;
    }
    position *= 1.0/w_sum;
    bump.Position = position;
}

static bump_t merge_bumps(const std::vector<bump_t> bumps) {
    bump_t bump = bumps[0];
    for(size_t i=1;i<bumps.size();i++) {
        bump_t b = bumps[i];
        for(size_t j=0;j<bump.Weights.size();j++) {
            bump.Weights[j] += b.Weights[j];
        }
    }
    // normalize
    double w_max = 0;
    size_t i_max = 0;
    for(size_t i=0;i<bump.Weights.size();i++) {
        bump.Weights[i] /= bumps.size();
        if(w_max<bump.Weights[i]) {
            i_max = i;
            w_max = bump.Weights[i];
        }
    }
    bump.MaxIndex = i_max;
    return bump;
}

static void split_cluster(const std::vector<crystal_t>& cluster,
                          std::vector< std::vector<crystal_t> >& clusters) {

    // make Voting based on relative distance or energy difference

    double totalClusterEnergy = 0;
    std::vector<unsigned> votes;
    votes.resize(cluster.size(), 0);
    // start searching at the second highest energy (i>0 case in next for loop)
    // since we know that the highest energy always has a vote
    votes[0]++;
    for(size_t i=0;i<cluster.size();i++) {
        totalClusterEnergy += cluster[i].Energy; // side calculation in this loop, but include i=0
        if(i==0)
            continue;

        // i>0 now...
        // for each crystal walk through cluster
        // according to energy gradient
        unsigned currPos = i;
        bool reachedMaxEnergy = false;
        double maxEnergy = 0;
        while(!reachedMaxEnergy) {
            // find neighbours intersection with actually hit clusters
            reachedMaxEnergy = true;
            const std::vector<unsigned>& neighbours = cluster[currPos].Element->Neighbours;
            for(size_t j=0;j<cluster.size();j++) {
                for(unsigned n=0;n<neighbours.size();n++) {
                    if(neighbours[n] != cluster[j].Element->Channel)
                        continue; // cluster element j not neighbour of element currPos, go to next
                    double energy = cluster[j].Energy;
                    if(maxEnergy < energy) {
                        maxEnergy = energy;
                        currPos = j;
                        reachedMaxEnergy = false;
                    }
                    break; // neighbour indices are unique, stop iterating
                }
            }
        }
        // currPos is now at max Energy
        votes[currPos]++;
    }

    // all crystals vote for highest energy
    // so this cluster should not be splitted,
    // just add it to the clusters
    if(votes[0] == cluster.size()) {
        clusters.push_back(cluster);
        return;
    }

    // find the bumps (crystals voted for)
    // and init the weights
    using bumps_t = std::list<bump_t>;
    bumps_t bumps;
    for(size_t i=0;i<votes.size();i++) {
        if(votes[i]==0)
            continue;
        // initialize the weights with the position of the crystal
        bump_t bump;
        bump.Position = cluster[i].Element->Position;
        bump.Weights.resize(cluster.size(), 0);
        calc_bump_weights(cluster, bump);
        bumps.push_back(bump);
    }

    // as long as we have overlapping bumps
    Bool_t haveOverlap = kFALSE;

    do {
        // converge the positions of the bumps
        unsigned iterations = 0;
        bumps_t stable_bumps;
        const double positionEpsilon = 0.01;
        while(!bumps.empty()) {
            for(bumps_t::iterator b=bumps.begin(); b != bumps.end(); ++b) {
                // calculate new bump position with current weights
                TVector3 oldPos = (*b).Position;
                update_bump_position(cluster, *b);
                double diff = (oldPos - (*b).Position).Mag();
                // check if position is stable
                if(diff>positionEpsilon) {
                    // no, then calc new weights with new position
                    calc_bump_weights(cluster, *b);
                    continue;
                }
                // yes, then save it and erase it from to-be-stabilized bumps
                stable_bumps.push_back(*b);
                b = bumps.erase(b);
            }
            // check max iterations
            iterations++;
            if(iterations<100)
                continue; // go on iterating...
            // ... or we iterated really long without convergence
            // discard the bumps which havn't converged
            bumps.clear();
        }

        // do we have any stable bumps?
        // Then just the use cluster as is
        if(stable_bumps.empty()) {
            clusters.push_back(cluster);
            return;
        }


        // stable_bumps are now identified, form clusters out of it
        // check if two bumps share the same crystal of highest energy
        // if they do, merge them

        typedef std::vector< std::vector< bump_t > > overlaps_t;
        overlaps_t overlaps(cluster.size()); // index of highest energy crystal -> corresponding stable bumps
        haveOverlap = kFALSE;
        for(bumps_t::iterator b=stable_bumps.begin(); b != stable_bumps.end(); ++b) {
            // remember the bump at its highest energy
            overlaps[b->MaxIndex].push_back(*b);
        }

        for(overlaps_t::iterator o=overlaps.begin(); o != overlaps.end(); ++o) {
            if(o->size()==0) {
                continue;
            }
            else if(o->size()==1) {
                bumps.push_back(o->at(0));
            }
            else { // size>1 more than one bump at index, then merge overlapping bumps
                haveOverlap = kTRUE;
                bumps.push_back(merge_bumps(*o));
            }
        }

    } while(haveOverlap);

    // bumps contain non-overlapping, stable bumps
    // try to build clusters out of it
    // we start with seeds at the position of the heighest weight in each bump,
    // and similarly to build_cluster iterate over the cluster's crystals


    // populate seeds and flags
    using bump_seeds_t = std::vector< std::vector<size_t> >;
    bump_seeds_t b_seeds; // for each bump, we track the seeds independently
    b_seeds.reserve(bumps.size());
    using state_t = std::vector< std::set<size_t> >;
    state_t state(cluster.size()); // at each crystal, we track the bump index
    for(bumps_t::iterator b=bumps.begin(); b != bumps.end(); ++b) {
        size_t i = b_seeds.size();
        state[b->MaxIndex].insert(i);
        // starting seed is just the max index
        std::vector<size_t> single;
        single.push_back(b->MaxIndex);
        b_seeds.push_back(single);
    }

    Bool_t noMoreSeeds = kFALSE;
    while(!noMoreSeeds) {
        bump_seeds_t b_next_seeds(bumps.size());
        state_t next_state = state;
        noMoreSeeds = kTRUE;
        for(size_t i=0; i<bumps.size(); i++) {
            // for each bump, do next neighbour iteration
            // so find intersection of neighbours of seeds with crystals inside the cluster
            std::vector<size_t> seeds = b_seeds[i];
            for(size_t j=0;j<cluster.size();j++) {
                // skip crystals in cluster which have already been visited/assigned
                if(state[j].size()>0)
                    continue;
                for(size_t s=0; s<seeds.size(); s++) {
                    crystal_t seed = cluster[seeds[s]];
                    for(size_t n=0;n<seed.Element->Neighbours.size();n++) {
                        if(seed.Element->Neighbours[n] != cluster[j].Element->Channel)
                            continue;
                        // for bump i, we found a next_seed, ...
                        b_next_seeds[i].push_back(j);
                        // ... and we assign it to this bump
                        next_state[j].insert(i);
                        // flag that we found more seeds
                        noMoreSeeds = kFALSE;
                        // neighbours is a list of unique items, we can stop searching
                        break;
                    }
                }
            }
        }

        // prepare for next iteration
        state = next_state;
        b_seeds = b_next_seeds;
    }

    // now, state tells us which crystals can be assigned directly to each bump
    // crystals are shared if they were claimed by more than one bump at the same neighbour iteration

    // first assign easy things and determine rough bump energy
    std::vector< std::vector<crystal_t> > bump_clusters(bumps.size());
    std::vector< double > bump_energies(bumps.size(), 0);
    for(size_t j=0;j<cluster.size();j++) {
        if(state[j].size()==1) {
            // crystal claimed by only one bump
            size_t i = *(state[j].begin());
            bump_clusters[i].push_back(cluster[j]);
            bump_energies[i] += cluster[j].Energy;
        }
    }

    // then calc weighted bump_positions for those preliminary bumps
    std::vector<TVector3> bump_positions(bumps.size(), TVector3(0,0,0));
    for(size_t i=0; i<bump_clusters.size(); i++) {
        std::vector<crystal_t> bump_cluster = bump_clusters[i];
        double w_sum = 0;
        for(size_t j=0;j<bump_cluster.size();j++) {
            double w = calc_energy_weight(bump_cluster[j].Energy, bump_energies[i]);
            bump_positions[i] += w * bump_cluster[j].Element->Position;
            w_sum += w;
        }
        bump_positions[i] *= 1.0/w_sum;
    }

    // finally we can share the energy of crystals claimed by more than one bump
    // we use bump_positions and bump_energies to do that
    for(size_t j=0;j<cluster.size();j++) {
        if(state[j].size()==1)
            continue;
        // size should never be zero, aka a crystal always belongs to at least one bump

        std::vector<double> pulls(bumps.size());
        double sum_pull = 0;
        for(std::set<size_t>::iterator b=state[j].begin(); b != state[j].end(); ++b) {
            TVector3 r = cluster[j].Element->Position - bump_positions[*b];
            double pull = bump_energies[*b] * TMath::Exp(-r.Mag()/cluster[j].Element->MoliereRadius);
            pulls[*b] = pull;
            sum_pull += pull;
        }

        for(std::set<size_t>::iterator b=state[j].begin(); b != state[j].end(); ++b) {
            crystal_t crys = cluster[j];
            crys.Energy *= pulls[*b]/sum_pull;
            bump_clusters[*b].push_back(crys);
        }
    }

    for(size_t i=0; i<bump_clusters.size(); i++) {
        std::vector<crystal_t> bump_cluster = bump_clusters[i];
        // always sort before adding to clusters
        sort(bump_cluster.begin(), bump_cluster.end());
        clusters.push_back(bump_cluster);
    }
}

static void build_cluster(std::list<crystal_t>& crystals,
                          std::vector<crystal_t>& cluster) {
    // first crystal has highest energy
    std::list<crystal_t>::iterator i = crystals.begin();

    // start with initial seed list
    std::vector<crystal_t> seeds;
    seeds.push_back(*i);

    // save i in the current cluster
    cluster.push_back(*i);
    // remove it from the candidates
    crystals.erase(i);

    while(seeds.size()>0) {
        // neighbours of all seeds are next seeds
        std::vector<crystal_t> next_seeds;

        for(std::vector<crystal_t>::iterator seed=seeds.begin(); seed != seeds.end(); seed++) {
            // find intersection of neighbours and seed
            for(std::list<crystal_t>::iterator j = crystals.begin() ; j != crystals.end() ; ) {
                bool foundNeighbour = false;
                for(size_t n=0;n<(*seed).Element->Neighbours.size();n++) {
                    if((*seed).Element->Neighbours[n] != (*j).Element->Channel)
                        continue;
                    next_seeds.push_back(*j);
                    cluster.push_back(*j);
                    j = crystals.erase(j);
                    foundNeighbour = true;
                    // neighbours is a list of unique items, we can stop searching
                    break;
                }
                // removal moves iterator already one forward
                if(!foundNeighbour)
                    ++j;
            }
        }
        // set new seeds, if any new found...
        seeds = next_seeds;
    }

    // sort it by energy
    sort(cluster.begin(), cluster.end());
}

static void do_clustering(
        std::list<crystal_t>& crystals,
        std::vector< std::vector<crystal_t> >& clusters
        ) {
    crystals.sort();

    while(crystals.size()>0) {
        std::vector<crystal_t> cluster;
        build_cluster(crystals, cluster); // already sorts "cluster" it by energy
        split_cluster(cluster, clusters);
    }
}


}}} // ant::reconstruct::clustering
