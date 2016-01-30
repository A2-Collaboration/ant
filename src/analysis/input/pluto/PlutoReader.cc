#include "PlutoReader.h"

#include <string>
#include <iostream>
#include <memory>

#include "utils/particle_tools.h"

#include "detail/PlutoWrapper.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"

#include "TTree.h"
#include "TClonesArray.h"

#include "expconfig/ExpConfig.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;
using namespace ant::analysis::data;
using namespace std;

PlutoReader::PlutoReader(const std::shared_ptr<WrapTFileInput>& rootfiles) :
    files(rootfiles)
{
    try {
        tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    }
    catch(ExpConfig::Exception e) {
        LOG(WARNING) << "Not generating MCTrue tagger hits since no tagger detector was found: " << e.what();
    }

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

    LOG(INFO) << "MCTrue input active, " << (tid_from_file ? "with" : " WITHOUT") << " TID match check";
}

PlutoReader::~PlutoReader() {}


const ParticleTypeDatabase::Type* PlutoReader::GetType(const PParticle* p) const {

    const ParticleTypeDatabase::Type* type= ParticleTypeDatabase::GetTypeOfPlutoID( p->ID() );

    if(!type) {

        type = ParticleTypeDatabase::AddTempPlutoType(
                    p->ID(),
                    "Pluto_"+to_string(p->ID()),
                    "Pluto_"+ string( pluto_database->GetParticleName(p->ID()) ),
                    pluto_database->GetParticleMass(p->ID())*1000.0,
                    pluto_database->GetParticleCharge(p->ID()) != 0
                );
        VLOG(7) << "Adding a ParticleType for pluto ID " << p->ID() << " to ParticleTypeDatabse";

        if(!type)
            /// \todo Change this so that it does not stop the whole process and can be caught for each particle
            throw std::out_of_range("Could not create dynamic mapping for Pluto Particle ID "+to_string(p->ID()));
    }

    return type;
}

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

void PlutoReader::CopyPluto(TEvent::Data& mctrue, const TID& eventid_rec)
{
    const Long64_t entries = PlutoMCTrue->GetEntries();
    std::vector<const PParticle*> PlutoParticles;
    PlutoParticles.reserve(size_t(entries));

    std::vector<size_t> dilepton_indices;

    for( Long64_t i=0; i < entries; ++i) {

        const PParticle* particle = dynamic_cast<const PParticle*>((*PlutoMCTrue)[i]);

        // ignore those weird dilepton particles
        if(particle->ID() == pluto_database->GetParticleID("dilepton"))
            dilepton_indices.push_back(i);
        PlutoParticles.push_back(particle);
    }

    using treenode_t = shared_ptr<Tree<TParticlePtr>>;

    vector<treenode_t> FlatTree;

    // convert pluto particles to ant particles and place in buffer list

    for( auto& PlutoParticle : PlutoParticles ) {

        // make an AntParticle out of it
        auto type = GetType(PlutoParticle);
        TLorentzVector lv = *PlutoParticle;
        lv *= 1000.0;   // convert to MeV
        auto AntParticle = make_shared<TParticle>(*type,lv);

        // Add particle to event storage
        if(PlutoParticle->GetDaughterIndex() == -1 ) { // final state
            mctrue.Particles.Add(AntParticle);
        }

        // Simulate some tagger hit
        if( AntParticle->Type() == ParticleTypeDatabase::BeamTarget) {
            const double energy = AntParticle->Ek();
            unsigned channel = 0;
            if(tagger && tagger->TryGetChannelFromPhoton(energy, channel)) {
                const double time = 0.0; /// @todo handle non-prompt hits?
                mctrue.Tagger.Hits.emplace_back(channel, energy, time);
            }
        }

        // save it in list with same order as in PlutoParticles
        FlatTree.emplace_back(Tree<TParticlePtr>::MakeNode(AntParticle));

    }

    // build the particle tree

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

            if( TreeNode->Get()->Type() == ParticleTypeDatabase::BeamTarget) {
                auto& headnode = mctrue.ParticleTree;
                if(headnode) {
                    LOG(WARNING) << "Found more than one BeamTarget in MCTrue";
                }
                headnode = TreeNode;
            }
            else
            {
                auto search_result = FindParticleByID(PlutoParticles, PlutoParticle->GetParentId());

                if(search_result.second) {
                    VLOG(7) << "Recovered missing pluto decay tree information.";
                    TreeNode->SetParent(FlatTree.at(search_result.first));
                } else {  // BeamProton is not supposed to have a parent
                    missing_decay_treeinfo = true;
                }

            }
        }
    }



    // for gun generated pluto things, there's no BeamTarget particle and thus no tree...
    if(mctrue.ParticleTree) {
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
        mctrue.ParticleTree->Sort(utils::ParticleTools::SortParticleByName);

    }

    // use eventID from file if available
    // also check if it matches with reconstructed TID
    if(tid_from_file) {
        mctrue.ID = *tid;
        if(!eventid_rec.IsInvalid() && eventid_rec != mctrue.ID) {
            throw Exception(std_ext::formatter()
                            << "TID mismatch: Reconstructed=" << eventid_rec
                            << " not equal to MCTrue=" << mctrue.ID);
        }
    }

    // calculate energy sum based on direction of particle
    auto& triggerinfos = mctrue.Trigger;
    double Esum = 0;
    for(const TParticlePtr& particle : mctrue.Particles.GetAll()) {
        if(geometry.DetectorFromAngles(*particle) & Detector_t::Type_t::CB) {
            Esum += particle->Ek();
        }
    }
    triggerinfos.CBEnergySum = Esum;

    /// @note multiplicity is only known on reconstructed

    // dump some information about the conversion

    if(missing_decay_treeinfo && mctrue.ParticleTree) {
        const auto& tid = !mctrue.ID.IsInvalid() ? mctrue.ID : eventid_rec;
        LOG(WARNING) << "Missing decay tree info for event " << tid;
        VLOG(5)      << "Dumping Pluto particles:\n" << PlutoTable(PlutoParticles);
        mctrue.ParticleTree = nullptr;
    }

    if(!dilepton_indices.empty()) {
        VLOG(5) << "Particle tree cleaned from dileptons: " << utils::ParticleTools::GetDecayString(mctrue.ParticleTree);
        VLOG(5) << "Dumping Pluto particles:\n" << PlutoTable(PlutoParticles);
    }
}


bool PlutoReader::ReadNextEvent(TEvent& event)
{
    if(!tree)
        return false;

    if(current_entry >= tree->GetEntries())
        return false;

    tree->GetEntry(current_entry);

    CopyPluto(*event.MCTrue, event.Reconstructed->ID);

    ++current_entry;

    return true;
}

double PlutoReader::PercentDone() const
{
    return double(current_entry) / double(tree->GetEntries());
}
