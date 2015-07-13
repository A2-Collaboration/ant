#include "GoatReader.h"

#include <string>
#include <iostream>
#include <memory>

#include "TTree.h"

#include "detail/PlutoWrapper.h"

#include "base/Logger.h"

using namespace ant;
using namespace ant::input;
using namespace std;

/**
 * @brief map goat apparatus numbers to apparatus_t enum values
 * in case unknown values show up: -> exception and do not sliently ignore
 */
detector_t IntToDetector_t( const int& a ) {
    detector_t d = detector_t::None;
    if(a & TrackInput::DETECTOR_NaI) {
        d |= detector_t::NaI;
    }
    if(a & TrackInput::DETECTOR_PID) {
        d |= detector_t::PID;
    }
    if(a & TrackInput::DETECTOR_MWPC) {
        d |= detector_t::MWPC;
    }
    if(a & TrackInput::DETECTOR_BaF2) {
        d |= detector_t::BaF2;
    }
    if(a & TrackInput::DETECTOR_PbWO4) {
        d |= detector_t::PbWO4;
    }
    if(a & TrackInput::DETECTOR_Veto) {
        d |= detector_t::Veto;
    }
    return d;
}

void GoatReader::CopyTagger(std::shared_ptr<Event> &event)
{
    for( Int_t i=0; i<tagger.GetNTagged(); ++i) {
        event->Reconstructed().TaggerHits().emplace_back(
                    TaggerHitPtr(new TaggerHit(
                                     tagger.GetTaggedChannel(i),
                                     tagger.GetTaggedEnergy(i),
                                     tagger.GetTaggedTime(i))
                                 ));
    }
}

void GoatReader::CopyTrigger(std::shared_ptr<Event> &event)
{
    TriggerInfo& ti = event->Reconstructed().TriggerInfos();

    ti.CBEenergySum() = trigger.GetEnergySum();
    ti.Multiplicity() = trigger.GetMultiplicity();

    for( int err=0; err < trigger.GetNErrors(); ++err) {
        ti.Errors().emplace_back(
                    DAQError(
                        trigger.GetErrorModuleID()[err],
                        trigger.GetErrorModuleIndex()[err],
                        trigger.GetErrorCode()[err]));
    }
}

void GoatReader::CopyDetectorHits(std::shared_ptr<Event>&)
{

}

/**
 * @brief map the cluster sizes from goat to unisgend ints
 * negative values mean no hit in the calorimeter
 * map those to 0
 */
clustersize_t GoatReader::MapClusterSize(const int& size) {
    return size < 0 ? 0 : size;
}

void GoatReader::CopyTracks(std::shared_ptr<Event> &event)
{
    for(Int_t i=0; i< tracks.GetNTracks(); ++i) {

        event->Reconstructed().Tracks().emplace_back(
                    TrackPtr( new Track(
                                  tracks.GetClusterEnergy(i),
                                  tracks.GetTheta(i),
                                  tracks.GetPhi(i),
                                  tracks.GetTime(i),
                                  MapClusterSize(tracks.GetClusterSize(i)),
                                  IntToDetector_t(tracks.GetDetectors(i)),
                                  tracks.GetVetoEnergy(i),
                                  tracks.GetMWPC0Energy(i),
                                  tracks.GetMWPC1Energy(i)
                                  )));
    }
}

const ParticleTypeDatabase::Type* GoatReader::GetType(const PParticle* p) const {

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
            // TODO: Change this so that it does not stop the whole process and can be caught for each particle
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

void GoatReader::CopyPluto(std::shared_ptr<Event> &event)
{
    const auto particles = pluto.Particles();

    ParticleList buffer;
    buffer.reserve(particles.size());

    // convert pluto particles to ant particles and place in buffer list
    for( auto& p : particles ) {

        auto type = GetType(p);
        TLorentzVector lv = *p;
        lv *= 1000.0;   // convert to MeV

        buffer.emplace_back(new Particle(*type,lv));

    }

    // loop over both lists (pluto and ant particles)
    for(size_t i=0; i<particles.size(); ++i) {

        const PParticle* plp = particles.at(i);
        ParticlePtr      alp = buffer.at(i);

        // Set up tree relations (parent/daughter)
        if(plp->GetParentIndex() >= 0) {

            if( size_t(plp->GetParentIndex()) < buffer.size()) {
                ParticlePtr parent = buffer.at(plp->GetParentIndex());
                alp->Partents().emplace_back(parent);
                parent->Daughters().emplace_back(alp);
            }

        } else {
            if(! (alp->Type() == ParticleTypeDatabase::BeamProton)) {

                auto sr = FindParticleByID(particles, plp->GetParentId());

                if(sr.second) {
                    VLOG(7) << "Recovered missing pluto decay tree inforamtion.";
                    ParticlePtr parent = buffer.at(sr.first);
                    alp->Partents().emplace_back(parent);
                    parent->Daughters().emplace_back(alp);
                } else {
                    LOG(WARNING) << "Missing decay tree info for pluto particle";  // BeamProton is not supposed to have a parent
                }
            }
        }

        // create a tagger hit corresponding to beam+parget particle
        if( alp->Type() == ParticleTypeDatabase::BeamProton) {
            const double energy = alp->E() - ParticleTypeDatabase::Proton.Mass();
            const int channel = 0; //TODO: Get tagger channel from energy -> Tagger cfg
            const double time = 0.0;
            event->MCTrue().TaggerHits().emplace_back( TaggerHitPtr( new TaggerHit(channel, energy, time) ) );
        }

        // Add particle to event storage
        if(plp->GetDaughterIndex() == -1 ) { // final state
            event->MCTrue().Particles().AddParticle(alp);
        } else { //intermediate
            event->MCTrue().Intermediates().AddParticle(alp);
        }
    }
    //TODO: CBEsum/Multiplicity into TriggerInfo
}

void GoatReader::CopyParticles(std::shared_ptr<Event> &event, ParticleInput &input_module, const ParticleTypeDatabase::Type &type)
{
    for(Int_t i=0; i < input_module.GetNParticles(); ++i) {

        const auto trackIndex = input_module.GetTrackIndex(i);
        if(trackIndex == -1) {
            cerr << "No Track for this particle!!" << endl;
        } else {
            const auto& track = event->Reconstructed().Tracks().at(trackIndex);

            event->Reconstructed().Particles().AddParticle(
                    std::make_shared<Particle>(type,track));
        }

    }
}

GoatReader::GoatReader(): pluto_database(makeStaticData())
{
}

void GoatReader::AddInputFile(const std::string &filename)
{
    files.OpenFile(filename);
}

class MyTreeRequestMgr: public TreeRequestManager {
protected:
    FileManager& m;
    TreeManager& t;

public:
    MyTreeRequestMgr(FileManager& _m, TreeManager& _t):
        m(_m), t(_t) {}

    TTree *GetTree(const std::string &name) {
        TTree* tree = nullptr;
        if( m.GetObject(name, tree) ) {
            t.AddTree(tree);
            VLOG(6) << "TTree " << name << " opened";
            return tree;
        } else
            VLOG(7) << "Could not find TTree " << name << " in any open GoAT file";
            return nullptr;
    }

};

//TODO: find a smart way to manage trees and modules:
//   if module does not init or gets removed-> remove also the tree from the list
//   two modules use same tree?
//   reset branch addresses ?

void GoatReader::Initialize()
{
    for(auto module = active_modules.begin(); module != active_modules.end(); ) {

        if( (*module)->SetupBranches( MyTreeRequestMgr(files, trees))) {
            module++;
        } else {
            module = active_modules.erase(module);
            VLOG(7) << "Not activating GoAT Input module";
        }
    }

}

Long64_t GoatReader::GetNEvents() const
{
    return trees.GetEntries();
}

bool GoatReader::hasData() const {
    return current_entry+1 < GetNEvents();
}

long long GoatReader::EventsRead() const
{
    return current_entry;
}

long long GoatReader::TotalEvents() const
{
    return GetNEvents();
}

std::shared_ptr<Event> GoatReader::ReadNextEvent()
{
    ++current_entry;
    trees.GetEntry(current_entry);

    active_modules.GetEntry();

    std::shared_ptr<Event> event = make_shared<Event>();

    CopyTrigger(event);
    CopyTagger(event);
    CopyTracks(event);
    CopyPluto(event);
    CopyDetectorHits(event);

    CopyParticles(event, photons, ParticleTypeDatabase::Photon);
    CopyParticles(event, protons, ParticleTypeDatabase::Proton);
    CopyParticles(event, pichagred, ParticleTypeDatabase::eCharged);
    CopyParticles(event, echarged, ParticleTypeDatabase::PiCharged);
    CopyParticles(event, neutrons, ParticleTypeDatabase::Neutron);


    return event;
}
