#include "GoatReader.h"

#include "detail/TreeManager.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/memory.h"

#include "TTree.h"

#include <string>
#include <iostream>
#include <memory>


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;
using namespace ant::analysis::data;
using namespace std;

/**
 * @brief map goat apparatus numbers to apparatus_t enum values
 * in case unknown values show up: -> exception and do not sliently ignore
 */
Detector_t::Any_t IntToDetector_t(const int& a) {
    auto d = Detector_t::Any_t::None;
    if(a & TrackInput::DETECTOR_NaI) {
        d |= Detector_t::Type_t::CB;
    }
    if(a & TrackInput::DETECTOR_PID) {
        d |= Detector_t::Type_t::PID;
    }
    if(a & TrackInput::DETECTOR_MWPC) {
        d |= Detector_t::Any_t::Tracker;
    }
    if(a & TrackInput::DETECTOR_BaF2) {
        d |= Detector_t::Type_t::TAPS;
    }
    if(a & TrackInput::DETECTOR_PbWO4) {
        d |= Detector_t::Type_t::TAPS;
    }
    if(a & TrackInput::DETECTOR_Veto) {
        d |= Detector_t::Type_t::TAPSVeto;
    }
    return d;
}

void GoatReader::CopyTagger(Event& event)
{
    for( Int_t i=0; i<tagger.GetNTagged(); ++i) {
        event.Reconstructed.TaggerHits.emplace_back(
                    tagger.GetTaggedChannel(i),
                    tagger.GetTaggedEnergy(i),
                    tagger.GetTaggedTime(i)
                    );
    }
}

void GoatReader::CopyTrigger(Event& event)
{
    Trigger_t& ti = event.Reconstructed.Trigger;

    ti.CBEnergySum = trigger.GetEnergySum();
    ti.ClusterMultiplicity = trigger.GetMultiplicity();

    for( int err=0; err < trigger.GetNErrors(); ++err) {
        ti.Errors.emplace_back(
                    DAQError(
                        trigger.GetErrorModuleID()[err],
                        trigger.GetErrorModuleIndex()[err],
                        trigger.GetErrorCode()[err]));
    }
}


/**
 * @brief map the cluster sizes from goat to unisgend ints
 * negative values mean no hit in the calorimeter
 * map those to 0
 */
clustersize_t GoatReader::MapClusterSize(const int& size) {
    return size < 0 ? 0 : size;
}

void GoatReader::CopyTracks(Event& event)
{
    for(Int_t i=0; i< tracks.GetNTracks(); ++i) {

        const auto det = IntToDetector_t(tracks.GetDetectors(i));

        Detector_t::Type_t d = Detector_t::Type_t::CB;

        if(det & Detector_t::Type_t::TAPS) {
            d = Detector_t::Type_t::TAPS;
        }

        event.Reconstructed.Candidates.emplace_back(
                    make_shared<TCandidate>(
                        det,
                        tracks.GetClusterEnergy(i),
                        tracks.GetTheta(i),
                        tracks.GetPhi(i),
                        tracks.GetTime(i),
                        MapClusterSize(tracks.GetClusterSize(i)),
                        tracks.GetVetoEnergy(i),
                        tracks.GetMWPC0Energy(i)+tracks.GetMWPC1Energy(i),
                        std::vector<TCluster>{ TCluster(TVector3(),tracks.GetClusterEnergy(i),tracks.GetTime(i),d,0)
                        }
                        )
                    );
    }
}


void GoatReader::CopyParticles(Event& event, ParticleInput& input_module,
                               const ParticleTypeDatabase::Type& type)
{
    for(Int_t i=0; i < input_module.GetNParticles(); ++i) {

        const auto trackIndex = input_module.GeTCandidateIndex(i);
        if(trackIndex == -1) {
            cerr << "No Track for this particle!!" << endl;
        } else {
            const auto& track = event.Reconstructed.Candidates.at(trackIndex);

            event.Reconstructed.Particles.AddParticle(
                    std::make_shared<Particle>(type,track));
        }

    }
}

class MyTreeRequestMgr: public TreeRequestManager {
protected:
    WrapTFileInput& m;
    TreeManager& t;

public:
    MyTreeRequestMgr(WrapTFileInput& _m, TreeManager& _t):
        m(_m), t(_t) {}

    TTree *GetTree(const std::string &name) {
        TTree* tree = nullptr;
        if( m.GetObject(name, tree) ) {
            t.AddTree(tree);
            VLOG(6) << "TTree " << name << " opened";
            return tree;
        } else
            VLOG(7) << "Could not find TTree " << name << " in any of the provided files";
            return nullptr;
    }

};

GoatReader::GoatReader(const std::shared_ptr<WrapTFileInput>& rootfiles):
    files(rootfiles),
    trees(std_ext::make_unique<TreeManager>())
{
    /// \todo find a smart way to manage trees and modules:
    //   if module does not init or gets removed-> remove also the tree from the list
    //   two modules use same tree?
    //   reset branch addresses ?

    for(auto module = active_modules.begin(); module != active_modules.end(); ) {

        if( (*module)->SetupBranches( MyTreeRequestMgr(*files, *trees))) {
            module++;
        } else {
            module = active_modules.erase(module);
            VLOG(7) << "Not activating GoAT Input module";
        }
    }
}

GoatReader::~GoatReader() {}

bool GoatReader::IsSource() { return trees->GetEntries()>0; }


bool GoatReader::ReadNextEvent(Event& event)
{
    if(current_entry>=trees->GetEntries())
        return false;

    trees->GetEntry(current_entry);

    active_modules.GetEntry();

    CopyTrigger(event);
    CopyTagger(event);
    CopyTracks(event);

    CopyParticles(event, photons, ParticleTypeDatabase::Photon);
    CopyParticles(event, protons, ParticleTypeDatabase::Proton);
    CopyParticles(event, pichagred, ParticleTypeDatabase::eCharged);
    CopyParticles(event, echarged, ParticleTypeDatabase::PiCharged);
    CopyParticles(event, neutrons, ParticleTypeDatabase::Neutron);

    ++current_entry;
    return true;
}

double GoatReader::PercentDone() const
{
    return double(current_entry)/double(trees->GetEntries());
}
