#include "ParticleID.h"

#include "tree/TParticle.h"

#include "base/std_ext/system.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TCutG.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;


ParticleID::~ParticleID() {}

TParticlePtr ParticleID::Process(const TCandidatePtr& cand) const
{
    auto type = Identify(cand);
    if(type !=nullptr) {
       return std::make_shared<TParticle>(*type, cand);
    }

    return nullptr;
}

std::unique_ptr<const ParticleID> ParticleID::default_particle_id = nullptr;

const ParticleID& ParticleID::GetDefault()
{
    return *default_particle_id;
}

void ParticleID::SetDefault(std::unique_ptr<const ParticleID> id)
{
    if(!default_particle_id)
        throw runtime_error("Default Particle ID has already been set and may not be changed anymore");
    default_particle_id = move(id);
}

SimpleParticleID::SimpleParticleID() {}

SimpleParticleID::~SimpleParticleID() {}

const ParticleTypeDatabase::Type* SimpleParticleID::Identify(const TCandidatePtr& cand) const
{
    if(cand->VetoEnergy>0.25)
        return addressof(ParticleTypeDatabase::Proton);
    else
        return addressof(ParticleTypeDatabase::Photon);
}

BasicParticleID::BasicParticleID() {}

BasicParticleID::~BasicParticleID()
{

}

bool TestCut(const std::shared_ptr<TCutG>& cut, const double& x, const double& y) {
    return (cut) && cut->IsInside(x,y);
}



const ParticleTypeDatabase::Type* BasicParticleID::Identify(const TCandidatePtr& cand) const
{
    const bool hadronic =    TestCut(tof,  cand->CaloEnergy, cand->Time)
                  || TestCut(size, cand->CaloEnergy, cand->ClusterSize);

    const bool hadronic_enabled = (tof) || (size);

    const bool charged = cand->VetoEnergy > 0.0;



    if(!charged) {

        if(hadronic) {
            return addressof(ParticleTypeDatabase::Neutron);
        } else {
            return addressof(ParticleTypeDatabase::Photon);
        }

    } else {  //charged

        if(
           (hadronic_enabled && hadronic)
           || (TestCut(dEE_proton, cand->CaloEnergy, cand->VetoEnergy))
           ) {
            return addressof(ParticleTypeDatabase::Proton);
        }

        if(
           TestCut(dEE_pion, cand->CaloEnergy, cand->VetoEnergy)
           ) {
            return addressof(ParticleTypeDatabase::PiCharged);
        }

        if(
           TestCut(dEE_electron, cand->CaloEnergy, cand->VetoEnergy)
           ) {
            return addressof(ParticleTypeDatabase::eCharged);
        }
    }

    return nullptr;
}




CBTAPSBasicParticleID::CBTAPSBasicParticleID(const string& pidcutsdir)
{
    try {
        WrapTFileInput cuts;
        VLOG(7) << "Looking for ParticleID cuts *.root in " << pidcutsdir;

        for(const auto& cutfile : std_ext::system::lsFiles(pidcutsdir,".root"))
        {
            try {
                cuts.OpenFile(cutfile);
            } catch (const std::runtime_error&) {
                LOG(WARNING) << "Could not open " << cutfile;
            }
        }
        LoadFrom(cuts);
    } catch (const std::runtime_error& e) {
        LOG(INFO) << "Failed to load cuts: " << e.what();
    }
}

CBTAPSBasicParticleID::~CBTAPSBasicParticleID()
{

}

const ParticleTypeDatabase::Type* CBTAPSBasicParticleID::Identify(const TCandidatePtr& cand) const
{
    if(cand->Detector & Detector_t::Any_t::CB_Apparatus) {
        return cb.Identify(cand);
    } else if(cand->Detector & Detector_t::Any_t::TAPS_Apparatus) {
        return taps.Identify(cand);
    }

    return nullptr;
}

void CBTAPSBasicParticleID::LoadFrom(WrapTFile& file)
{

        cb.dEE_proton       = file.GetSharedClone<TCutG>("cb_dEE_proton");
        cb.dEE_pion         = file.GetSharedClone<TCutG>("cb_dEE_pion");
        cb.dEE_electron     = file.GetSharedClone<TCutG>("cb_dEE_electron");
        cb.tof              = file.GetSharedClone<TCutG>("cb_ToF");
        cb.size             = file.GetSharedClone<TCutG>("cb_CluserSize");

        taps.dEE_proton     = file.GetSharedClone<TCutG>("taps_dEE_proton");
        taps.dEE_pion       = file.GetSharedClone<TCutG>("taps_dEE_pion");
        taps.dEE_electron   = file.GetSharedClone<TCutG>("taps_dEE_electron");
        taps.tof            = file.GetSharedClone<TCutG>("taps_ToF");
        taps.size           = file.GetSharedClone<TCutG>("taps_CluserSize");
}


