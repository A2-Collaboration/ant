#include "PlutoReader.h"

#include <string>
#include <iostream>
#include <memory>

#include "utils/particle_tools.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "detail/PlutoWrapper.h"

#include "expconfig/ExpConfig.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"
#include "base/std_ext/vector.h"

#include "TTree.h"
#include "TClonesArray.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

PlutoReader::PlutoReader(const std::shared_ptr<WrapTFileInput>& rootfiles) :
    files(rootfiles)
{
    /// \note the Pluto tree is really named just "data"
    if(!files->GetObject("data", tree))
        return;

    VLOG(5) << "Found Pluto Data Tree";

    const auto res = tree->SetBranchAddress("Particles",         &PlutoMCTrue);
    if(res != TTree::kMatch) LOG(ERROR) << "Could not access branch in PLuto tree";

    pluto_database = makeStaticData();

    TTree* tid_tree = nullptr;
    if(files->GetObject("data_tid", tid_tree)) {

        if(tid_tree->GetEntries() != tree->GetEntries()) {
            throw Exception("Pluto Tree / TID Tree size mismatch:");
        }

        tree->AddFriend(tid_tree);

        const auto res = tree->SetBranchAddress("tid", &tid);

        if(res==TTree::kMatch)
            tid_from_file = true;
    }
    else {
        /// \todo think of some better timestamp?
        tid = new TID(static_cast<std::uint32_t>(std::time(nullptr)),
                      0, // start with 0 as lower ID
                      std::list<TID::Flags_t>{TID::Flags_t::MC, TID::Flags_t::AdHoc} // mark as MC
                      );
        tid_from_file = false;
    }

    LOG(INFO) << "MCTrue input active" << (tid_from_file ? ", with TID match check" : "");
    LOG_IF(!tid_from_file, WARNING) << "No TID match check enabled";

    try {
        tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    }
    catch(ExpConfig::Exception e) {
        LOG(WARNING) << "Not generating MCTrue tagger hits since no tagger detector was found: " << e.what();
    }
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

void PlutoReader::CopyPluto(TEventData& mctrue)
{
    const Long64_t entries = PlutoMCTrue->GetEntries();
    std::vector<const PParticle*> PlutoParticles;
    PlutoParticles.reserve(size_t(entries));

    std::vector<size_t> dilepton_indices;

    for( Long64_t i=0; i < entries; ++i) {

        const PParticle* particle = dynamic_cast<const PParticle*>((*PlutoMCTrue)[i]);

        // remember positions of those weird dilepton/dimuon particles
        if(particle->ID() == pluto_database->GetParticleID("dilepton") ||
           particle->ID() == pluto_database->GetParticleID("dimuon"))
            dilepton_indices.push_back(i);
        PlutoParticles.push_back(particle);
    }

    using treenode_t = shared_ptr<Tree<TParticlePtr>>;

    vector<treenode_t> FlatTree;
    vector<shared_ptr<TParticle>> finalstate_particles;
    // convert pluto particles to ant particles and place in buffer list
    for(size_t i=0; i<PlutoParticles.size(); ++i) {

        auto& PlutoParticle = PlutoParticles[i];

        // find pluto type in database
        auto type = ParticleTypeDatabase::GetTypeOfPlutoID( PlutoParticle->ID() );

        // note that type might by nullptr (in particular for those dileptons...)
        // then just add some "empty" tree node
        if(!type) {
            if(!std_ext::contains(dilepton_indices, i))
                // check $PLUTOSYS/src/PStdData.cc what to do, you may add it to the ant Database if you wish
                throw Exception(std_ext::formatter() << "Unknown pluto particle found: ID="
                                << PlutoParticle->ID());

            FlatTree.emplace_back(Tree<TParticlePtr>::MakeNode(nullptr));
            continue;
        }

        // make an AntParticle out of it
        LorentzVec lv = *PlutoParticle;
        lv *= 1000.0;   // convert to MeV
        auto AntParticle = make_shared<TParticle>(*type,lv);

        // Consider final state particle
        if(PlutoParticle->GetDaughterIndex() == -1 ) { // final state
            finalstate_particles.emplace_back(AntParticle);
        }

        // Simulate some tagger hit
        if( AntParticle->Type() == ParticleTypeDatabase::BeamTarget) {
            const double energy = AntParticle->Ek();
            unsigned channel = 0;
            if(tagger && tagger->TryGetChannelFromPhoton(energy, channel)) {
                const double time = 0.0; /// @todo handle non-prompt hits?
                mctrue.TaggerHits.emplace_back(channel, energy, time);
            }
        }

        // save it in list with same order as in PlutoParticles
        FlatTree.emplace_back(Tree<TParticlePtr>::MakeNode(AntParticle));

    }

    // build the particle tree
    // for gun generated pluto things, the tree is just that one particle
    if(FlatTree.size() == 1) {
        mctrue.ParticleTree = FlatTree.front();
    }
    else {
        // loop over both lists (pluto and and the flat tree)
        bool missing_decay_treeinfo = false;
        for(size_t i=0; i<PlutoParticles.size(); ++i) {

            const PParticle* PlutoParticle = PlutoParticles.at(i);
            shared_ptr<Tree<TParticlePtr>>&  TreeNode = FlatTree.at(i);

            Int_t parent_index = PlutoParticle->GetParentIndex();

            // Set up tree relations (parent/daughter)
            if(parent_index >= 0 && size_t(parent_index) < FlatTree.size()) {

                TreeNode->SetParent(FlatTree.at(parent_index));

            } else {

                if(TreeNode->Get() &&
                   TreeNode->Get()->Type() == ParticleTypeDatabase::BeamTarget) {
                    auto& headnode = mctrue.ParticleTree;
                    if(headnode) {
                        LOG(WARNING) << "Found more than one BeamTarget in MCTrue, that's weird.";
                    }
                    headnode = TreeNode;
                }
                else
                {
                    auto search_result = FindParticleByID(PlutoParticles, PlutoParticle->GetParentId());

                    if(search_result.second) {
                        VLOG(7) << "Recovered missing pluto decay tree information.";
                        TreeNode->SetParent(FlatTree.at(search_result.first));
                    } else {
                        // recovery failed
                        missing_decay_treeinfo = true;
                    }

                }
            }
        }

        if(mctrue.ParticleTree) {
            if(missing_decay_treeinfo) {
                LOG(WARNING) << "Missing decay tree info for event " << mctrue.ID;
                VLOG(5)      << "Dumping Pluto particles:\n" << PlutoTable(PlutoParticles);
                mctrue.ParticleTree = nullptr;
            }
            else {
                // remove articifial Pluto_dilepton particles,
                // this assumes that a dilepton never
                // is a parent of a dilepton
                for(auto i : dilepton_indices) {
                    TParticleTree_t& dilepton = FlatTree[i];
                    while(!dilepton->Daughters().empty()) {
                        auto& d =  dilepton->Daughters().front();
                        d->SetParent(dilepton->GetParent());
                    }
                    dilepton->Unlink();
                }

                if(!dilepton_indices.empty()) {
                    VLOG(5) << "Particle tree cleaned from dileptons: " << utils::ParticleTools::GetDecayString(mctrue.ParticleTree);
                    VLOG(5) << "Dumping Pluto particles:\n" << PlutoTable(PlutoParticles);
                }

                mctrue.ParticleTree->Sort(utils::ParticleTools::SortParticleByName);
            }
        }
    }


    // calculate energy sum based on direction of particle
    auto& triggerinfos = mctrue.Trigger;
    triggerinfos.ClusterMultiplicity = 0;
    triggerinfos.CBTiming = 0;
    double Esum = 0;
    for(const auto& particle : finalstate_particles) {
        if(geometry.DetectorFromAngles(*particle) & Detector_t::Type_t::CB) {
            Esum += particle->Ek();
            // expected MCTrue multiplicity
            triggerinfos.ClusterMultiplicity++;
        }
    }
    triggerinfos.CBEnergySum = Esum;
}


bool PlutoReader::ReadNextEvent(event_t& event)
{
    if(!tree)
        return false;

    if(current_entry >= tree->GetEntries())
        return false;

    tree->GetEntry(current_entry);

    // use eventID from file if available
    // also check if it matches with reconstructed TID
    if(tid_from_file) {
        if(event.HasReconstructed() &&
           event.Reconstructed().ID != *tid) {
            throw Exception(std_ext::formatter()
                            << "TID mismatch: Reconstructed=" << event.Reconstructed().ID
                            << " not equal to MCTrue=" << *tid);
        }
    }

    // ensure MCTrue branch is there, potentially add TID if invalid so far
    if(!event.HasMCTrue()) {
        event.MakeMCTrue(*tid);
    }
    else if(event.MCTrue().ID.IsInvalid()) {
        event.MCTrue().ID = *tid;
    }

    CopyPluto(event.MCTrue());

    ++current_entry;

    if(!tid_from_file)
        ++(*tid);

    return true;
}

double PlutoReader::PercentDone() const
{
    return double(current_entry) / double(tree->GetEntries());
}
