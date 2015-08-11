#include "PlutoReader.h"

#include <string>
#include <iostream>
#include <memory>



#include "detail/PlutoWrapper.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext.h"

#include "TTree.h"
#include "TClonesArray.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;
using namespace ant::analysis::data;
using namespace std;

PlutoReader::PlutoReader(const std::shared_ptr<WrapTFileInput>& rootfiles):
    files(rootfiles)
{
    /// \note the Pluto tree is really named just "data"
    if(!files->GetObject("data", tree))
        return;

    tree->SetBranchAddress("Particles",         &PlutoMCTrue);
    tree->SetBranchAddress("plutoID",           &plutoID);
    tree->SetBranchAddress("plutoRandomID", 	&plutoRandomID);

    pluto_database = makeStaticData();
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
    particle->Parents().emplace_back(parent);
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
                    LOG(WARNING) << "Missing decay tree info for pluto particle";
                }

                VLOG(9) << "\n" << PlutoTable(PlutoParticles);
            }
        }

        // create a tagger hit corresponding to beam+parget particle
        if( AntParticle->Type() == ParticleTypeDatabase::BeamProton) {
            const double energy = AntParticle->E() - ParticleTypeDatabase::Proton.Mass();
            const int channel = 0; //TODO: Get tagger channel from energy -> Tagger cfg
            const double time = 0.0;
            event.MCTrue().TaggerHits().emplace_back( TaggerHitPtr( new TaggerHit(channel, energy, time) ) );
        }

        // Add particle to event storage
        if(PlutoParticle->GetDaughterIndex() == -1 ) { // final state
            event.MCTrue().Particles().AddParticle(AntParticle);
        } else { //intermediate
            event.MCTrue().Intermediates().AddParticle(AntParticle);
        }
    }

    /// \todo CBEsum/Multiplicity into TriggerInfo
}


bool PlutoReader::ReadNextEvent(Event& event, TSlowControl&)
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
