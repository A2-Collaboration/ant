#include "PlutoReader.h"

#include <string>
#include <iostream>
#include <memory>

#include "utils/ParticleTools.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "detail/PlutoWrapper.h"

#include "expconfig/ExpConfig.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"
#include "base/std_ext/container.h"

#include "TTree.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

PlutoReader::PlutoReader(const std::shared_ptr<WrapTFileInput>& rootfiles) :
    files(rootfiles)
{
    /// \note the Pluto tree is really named just "data"
    if(!files->GetObject("data", plutoTree.Tree))
        return;

    VLOG(5) << "Found Pluto 'data' tree";

    plutoTree.LinkBranches();

    if(files->GetObject("data_tid", tidTree.Tree)) {
        if(tidTree.Tree->GetEntries() != plutoTree.Tree->GetEntries()) {
            throw Exception("Pluto Tree / TID Tree size mismatch:");
        }
        tidTree.LinkBranches();
    }
    else {
        // think of some better timestamp here?
        tidTree.tid = TID(static_cast<std::uint32_t>(std::time(nullptr)),
                          0, // start with 0 as lower ID
                          std::list<TID::Flags_t>{TID::Flags_t::MC, TID::Flags_t::AdHoc} // mark as MC
                          );
    }

    LOG(INFO) << "MCTrue input active" << (tidTree ? ", with TID match check" : "") << ", entries=" << plutoTree.Tree->GetEntries();
    LOG_IF(!tidTree, WARNING) << "No TID match check enabled";

    try {
        tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    }
    catch(ExpConfig::Exception& e) {
        LOG(WARNING) << "Not generating MCTrue tagger hits since no tagger detector was found: " << e.what();
    }

    pluto_database = makeStaticData();
}

PlutoReader::~PlutoReader() {}

/**
 * @brief Find a PParticle in a vector by its ID. ID has to be unique in the vector.
 * @param particles vector to search
 * @param ID Id to look for
 * @return pair<size_t,bool>, first=Index of found particle, bool: true if found and unique, false else
 */
std::pair<std::size_t,bool> FindParticleByID(const std::vector<const PParticle*> particles, const int ID) {

    std::pair<std::size_t,bool> result = {0, false};

    for(std::size_t i = 0; i<particles.size(); ++i) {
        if(particles.at(i)->ID() == ID) {
            if(!result.second) {
                result.first = i;
                result.second = true;
            } else {
                result.second = false;
                return result;
            }
        }
    }

    return result;

}

std::string PlutoTable(const std::vector<const PParticle*> particles) {
    stringstream s;
    int i=0;
    s << "Index\tParticleName\tParentIndex\tNameOfParent\t|\tDaughterIndex\n";
    for(auto& p : particles ) {
        const PParticle* pl = p;
        s << i++ << "\t" << GetParticleName(p) << "\t" << p->GetParentIndex() << "\t" << GetParticleName(pl->GetParentId()) << "\t|\t" <<  p->GetDaughterIndex() << "\t" <<  endl;
    }

    return s.str();
}

using PlutoParticles_t = std::vector<const PParticle*>;

bool BuildParticleGunTree(
        TParticleTree_t& mctrueTree,
        const PlutoParticles_t& plutoParticles,
        const vector<TParticleTree_t>& flatTree)
{
    // any particle with parent or daughter
    // means it's not a gun
    for(auto& plutoParticle : plutoParticles) {
        if(plutoParticle->GetParentIndex()>=0)
            return false;
        if(plutoParticle->GetDaughterIndex()>=0)
            return false;
    }

    // BeamTarget type is also not allowed in guns
    for(auto& antParticleTree : flatTree) {
        if(antParticleTree->Get()->Type() == ParticleTypeDatabase::BeamTarget)
            return false;
    }

    // build a tree with all particles as leaves and
    // "pseudo" GunParticle as headnode

    mctrueTree = Tree<TParticlePtr>::MakeNode(
                     make_shared<TParticle>(
                         ParticleTypeDatabase::ParticleGun,
                         // make sure really uses that as some "true" particle...
                         LorentzVec({std_ext::NaN,std_ext::NaN,std_ext::NaN}, std_ext::NaN)));

    for(auto& treeNode : flatTree) {
        treeNode->SetParent(mctrueTree);
    }

    return true;
}


bool BuildDecayTree(
        TParticleTree_t& mctrueTree,
        const PlutoParticles_t& plutoParticles,
        const vector<TParticleTree_t>& flatTree,
        const std::vector<size_t>& dileptonIndices)
{
    // loop over both lists (pluto and the flat tree)
    // in parallel
    bool missing_decay_treeinfo = false;
    for(size_t i=0; i<plutoParticles.size(); ++i) {

        auto plutoParticle = plutoParticles.at(i);
        auto& treeNode = flatTree.at(i);

        auto parent_index = plutoParticle->GetParentIndex();

        // Set up tree relations (parent/daughter)
        if(parent_index >= 0 && size_t(parent_index) < flatTree.size()) {

            treeNode->SetParent(flatTree.at(parent_index));

        } else {

            // remember the TreeNode might be nullptr (see above)
            if(treeNode->Get()) {
                if(mctrueTree) {
                    // Found more than one possible headnode in MCTrue
                    mctrueTree = nullptr;
                    return false;
                }
                mctrueTree = treeNode;
            }
            else
            {
                auto search_result = FindParticleByID(plutoParticles, plutoParticle->GetParentId());

                if(search_result.second) {
                    VLOG(7) << "Recovered missing pluto decay tree information.";
                    treeNode->SetParent(flatTree.at(search_result.first));
                } else {
                    // recovery failed
                    missing_decay_treeinfo = true;
                }

            }
        }
    }

    if(!mctrueTree)
        return false;

    if(missing_decay_treeinfo) {
        mctrueTree = nullptr;
        return false;
    }

    // remove artificial Pluto_dilepton particles,
    // this assumes that a dilepton never
    // is a parent of a dilepton
    for(auto i : dileptonIndices) {
        auto& dilepton = flatTree[i];
        while(!dilepton->Daughters().empty()) {
            auto& d =  dilepton->Daughters().front();
            d->SetParent(dilepton->GetParent());
        }
        dilepton->Unlink();
    }

    mctrueTree->Sort(utils::ParticleTools::SortParticleByName);

    return true;
}

void PlutoReader::CopyPluto(TEventData& mctrue)
{
    const auto nParticles = plutoTree.Particles().GetEntries();
    PlutoParticles_t plutoParticles;
    plutoParticles.reserve(size_t(nParticles));

    std::vector<size_t> dileptonIndices;

    for(auto i=0;i<nParticles;++i) {
        auto particle = dynamic_cast<const PParticle*>(plutoTree.Particles()[i]);
        // remember positions of those weird dilepton/dimuon particles
        if(particle->ID() == pluto_database->GetParticleID("dilepton") ||
           particle->ID() == pluto_database->GetParticleID("dimuon"))
            dileptonIndices.push_back(i);
        plutoParticles.push_back(particle);
    }


    vector<TParticleTree_t> flatTree;
    TParticleList finalstateParticles;
    // convert pluto particles to ant particles and place in buffer list
    for(size_t i=0; i<plutoParticles.size(); ++i) {

        auto& plutoParticle = plutoParticles[i];

        // find pluto type in database
        auto type = ParticleTypeDatabase::GetTypeFromPlutoID( plutoParticle->ID() );

        // note that type might be nullptr (in particular for those dileptons...)
        // then just add some "empty" tree node
        if(!type) {
            if(!std_ext::contains(dileptonIndices, i))
                // check $PLUTOSYS/src/PStdData.cc what to do, you may add it to the ant Database if you wish
                throw Exception(std_ext::formatter() << "Unknown pluto particle found: ID="
                                << plutoParticle->ID());

            flatTree.emplace_back(Tree<TParticlePtr>::MakeNode(nullptr));
            continue;
        }

        // make an AntParticle out of it
        LorentzVec lv = *plutoParticle;
        lv *= 1000.0;   // convert to MeV
        auto antParticle = make_shared<TParticle>(*type,lv);

        // Consider final state particle
        if(plutoParticle->GetDaughterIndex() == -1 ) { // final state
            finalstateParticles.emplace_back(antParticle);
        }

        // Simulate some tagger hit
        if( antParticle->Type() == ParticleTypeDatabase::BeamTarget) {
            const double energy = antParticle->Ek();
            unsigned channel = 0;
            if(tagger && tagger->TryGetChannelFromPhoton(energy, channel)) {
                const double time = 0.0; /// @todo handle non-prompt hits?
                mctrue.TaggerHits.emplace_back(channel, energy, time);
            }
        }

        // save it in list with same order as in PlutoParticles
        flatTree.emplace_back(Tree<TParticlePtr>::MakeNode(antParticle));

    }

    // try building the particle tree
    // first check if it's a particle gun event
    // then try building the usual decay tree
    if(!BuildParticleGunTree(mctrue.ParticleTree, plutoParticles, flatTree)) {
        if(!BuildDecayTree(mctrue.ParticleTree, plutoParticles, flatTree, dileptonIndices)) {
            LOG_N_TIMES(10, WARNING) << "Missing decay tree info for event " << mctrue.ID
                                     << " (max 10 times reported)";
            VLOG(5)      << "Dumping Pluto particles:\n" << PlutoTable(plutoParticles);
        }
    }

    // some diagnostics
    if(!dileptonIndices.empty()) {
        VLOG(5) << "Particle tree cleaned from dileptons: " << utils::ParticleTools::GetDecayString(mctrue.ParticleTree);
        VLOG(5) << "Dumping Pluto particles:\n" << PlutoTable(plutoParticles);
    }

    // calculate energy sum based on direction of particle
    auto& triggerInfos = mctrue.Trigger;
    triggerInfos.ClusterMultiplicity = 0;
    triggerInfos.CBTiming = 0;
    double Esum = 0;
    for(const auto& particle : finalstateParticles) {
        if(geometry.DetectorFromAngles(*particle) & Detector_t::Type_t::CB) {
            Esum += particle->Ek();
            // expected MCTrue multiplicity
            triggerInfos.ClusterMultiplicity++;
        }
    }
    triggerInfos.CBEnergySum = Esum;
}


bool PlutoReader::ReadNextEvent(event_t& event)
{
    if(!plutoTree)
        return false;

    if(current_entry >= plutoTree.Tree->GetEntries())
        return false;

    plutoTree.Tree->GetEntry(current_entry);

    // use eventID from file if available
    // also check if it matches with reconstructed TID
    if(tidTree) {
        tidTree.Tree->GetEntry(current_entry);
        if(event.HasReconstructed() &&
           event.Reconstructed().ID != tidTree.tid) {
            throw Exception(std_ext::formatter()
                            << "TID mismatch: Reconstructed=" << event.Reconstructed().ID
                            << " not equal to MCTrue=" << tidTree.tid);
        }
    }

    // ensure MCTrue branch is there, potentially add TID if invalid so far
    if(!event.HasMCTrue()) {
        event.MakeMCTrue(tidTree.tid);
    }
    else if(event.MCTrue().ID.IsInvalid()) {
        event.MCTrue().ID = tidTree.tid;
    }

    CopyPluto(event.MCTrue());

    ++current_entry;

    if(!tidTree)
        ++tidTree.tid();

    return true;
}

double PlutoReader::PercentDone() const
{
    return double(current_entry) / double(plutoTree.Tree->GetEntries());
}
