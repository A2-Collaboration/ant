#include "PlutoReader.h"

#include <string>
#include <iostream>
#include <memory>



#include "detail/PlutoWrapper.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"

#include "TTree.h"
#include "TClonesArray.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;
using namespace ant::analysis::data;
using namespace std;

PlutoReader::PlutoReader(const std::shared_ptr<WrapTFileInput>& rootfiles, const std::shared_ptr<TaggerDetector_t> tagger_) :
    tagger(tagger_),
    files(rootfiles)
{
    if(!tagger) {
        LOG(WARNING) << "Not generating MCTrue tagger hits since no tagger detector was provided";
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
            throw Exception("Pluto Tree / TID Tree size missmatch:");
        }

        tree->AddFriend(tid_tree);

        const auto res = tree->SetBranchAddress("tid", &tid);

        if(res==TTree::kMatch)
            tid_from_file = true;
    }

    LOG(INFO) << "Pluto input active";
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

void SetParticleRelations(ParticleList list, ParticlePtr& particle, size_t parent_index) {
    ParticlePtr& parent = list.at(parent_index);
    particle->Parent() = parent;
    parent->Daughters().emplace_back(particle);
}

std::string PlutoTable(const std::vector<const PParticle*> particles) {
    stringstream s;
    int i=0;
    for(auto& p : particles ) {
        const PParticle* pl = p;
        s << i++ << "\t" << GetParticleName(p) << "\t" << p->GetParentIndex() << "\t" << GetParticleName(pl->GetParentId()) << "\t|\t" <<  p->GetDaughterIndex() << "\t" <<  endl;
    }

    return s.str();
}

void PlutoReader::CopyPluto(Event& event)
{
    const Long64_t entries = PlutoMCTrue->GetEntries();
    PlutoParticles.resize(0);
    PlutoParticles.reserve(size_t(entries));

    for( Long64_t i=0; i < entries; ++i) {

        const PParticle* const particle = dynamic_cast<const PParticle*>((*PlutoMCTrue)[i]);

        if(particle) {
            PlutoParticles.push_back(particle);
        }
    }

    ParticleList AntParticles;
    AntParticles.reserve(PlutoParticles.size());

    // convert pluto particles to ant particles and place in buffer list
    for( auto& p : PlutoParticles ) {

        auto type = GetType(p);
        TLorentzVector lv = *p;
        lv *= 1000.0;   // convert to MeV

        AntParticles.emplace_back(new Particle(*type,lv));

    }

    // loop over both lists (pluto and ant particles)
    for(size_t i=0; i<PlutoParticles.size(); ++i) {

        const PParticle* PlutoParticle = PlutoParticles.at(i);
        ParticlePtr      AntParticle = AntParticles.at(i);

        // Set up tree relations (parent/daughter)
        if(PlutoParticle->GetParentIndex() >= 0) {

            if( size_t(PlutoParticle->GetParentIndex()) < AntParticles.size()) {
                SetParticleRelations(AntParticles, AntParticle, PlutoParticle->GetParentIndex());
            }

        } else {
            if( ! (AntParticle->Type() == ParticleTypeDatabase::BeamProton)) {

                auto search_result = FindParticleByID(PlutoParticles, PlutoParticle->GetParentId());

                if(search_result.second) {
                    VLOG(7) << "Recovered missing pluto decay tree inforamtion.";
                    SetParticleRelations(AntParticles, AntParticle, search_result.first);

                } else {  // BeamProton is not supposed to have a parent
                    VLOG(7) << "Missing decay tree info for pluto particle";
                }

                VLOG(9) << "\n" << PlutoTable(PlutoParticles);
            }
        }

        // create a tagger hit corresponding to beam+parget particle
        if( AntParticle->Type() == ParticleTypeDatabase::BeamProton) {
            const double energy = AntParticle->E() - ParticleTypeDatabase::Proton.Mass();
            unsigned channel = 0;
            if(tagger && tagger->TryGetChannelFromPhoton(energy, channel)) {
                const double time = 0.0; /// @todo handle non-prompt hits?
                event.MCTrue().TaggerHits().emplace_back(make_shared<TaggerHit>(channel, energy, time));
            }
        }

        // Add particle to event storage
        if(PlutoParticle->GetDaughterIndex() == -1 ) { // final state
            event.MCTrue().Particles().AddParticle(AntParticle);
        } else { //intermediate
            event.MCTrue().Intermediates().AddParticle(AntParticle);
        }
    }

    auto& triggerinfos = event.MCTrue().TriggerInfos();

    // use eventID from file if available
    // also check if it matches with reconstructed TID
    if(tid_from_file) {
        triggerinfos.EventID() = *tid;

        const auto& eventid_rec = event.Reconstructed().TriggerInfos().EventID();
        if(!eventid_rec.IsInvalid() && eventid_rec != triggerinfos.EventID()) {
            throw Exception(std_ext::formatter()
                            << "TID mismatch: Reconstructed=" << eventid_rec
                            << " not equal to MCTrue=" << event.MCTrue().TriggerInfos().EventID());
        }
    }

    /// @note multiplicity is only known on reconstructed

}


bool PlutoReader::ReadNextEvent(Event& event)
{
    if(!tree)
        return false;

    if(current_entry >= tree->GetEntries())
        return false;

    tree->GetEntry(current_entry);

    CopyPluto(event);

    ++current_entry;

    return true;
}
