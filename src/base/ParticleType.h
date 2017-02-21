#pragma once

#include "base/types.h"
#include "base/interval.h"

#include <string>
#include <ostream>
#include <map>
#include <iterator>
#include <vector>

namespace ant {

class ParticleTypeDatabase {

    friend struct TParticle;

public:

    class Type {

        friend class ParticleTypeDatabase;
        friend struct TParticle;

    protected:
        unsigned UID; // used by TParticle for serialization
        static unsigned NextUID;
        std::string name;
        std::string print_name;
        mev_t mass;
        bool charged;
        const Type* sametype;

        Type(const std::string& _name, const std::string& _print_name, const mev_t& _mass, const bool& _charged, const Type* _sametype=nullptr );

    public:

        const std::string& PrintName()  const { return print_name; }
        const std::string& Name()       const { return name; }
        mev_t Mass()                    const { return mass; }
        bool  Charged()                 const { return charged; }
        interval<mev_t> GetWindow(const mev_t width) const { return {mass-width/2,mass+width/2}; }
        mev_t PhotoproductionThresh(const Type& target = ParticleTypeDatabase::Proton)   const { return CalculatePhotoproductionThreshold(mass, target); }

        bool operator== (const Type& rhs) const {
            if( this == &rhs )
                return true;
            if( sametype == &rhs )
                return true;
            if( rhs.sametype == this )
                return true;
            return false;
        }

        bool operator!= (const Type& rhs) const {
            return !this->operator==(rhs);
        }

        Type(const Type&) = delete;
        Type& operator=(const Type&) = delete;
    };

    using PlutoIDMap_t = std::map<index_t, const Type*>;
    using TypeList_t = std::vector<const Type*>;

protected:

    using  Particles_t = std::map<unsigned, const Type&>;
    static Particles_t types;

    static PlutoIDMap_t pluto_id_map;

    static const TypeList_t detectables;
    static const TypeList_t mc_finalstate;
    static const TypeList_t neutral_mesons;

public:

    ParticleTypeDatabase() {}

    virtual ~ParticleTypeDatabase() {}

    static void Print();

    static const Type Proton;
    static const Type Neutron;
    static const Type Nucleon;
    static const Type SigmaPlus;
    static const Type Photon;
    static const Type Pi0;
    static const Type PiPlus;
    static const Type PiMinus;
    static const Type PiCharged;

    static const Type K0s;

    static const Type ePlus;
    static const Type eMinus;
    static const Type eCharged;

    static const Type MuPlus;
    static const Type MuMinus;
    static const Type MuCharged;

    static const Type Eta;
    static const Type Omega;
    static const Type EtaPrime;
    static const Type Rho;

    static const Type BeamTarget;
    static const Type BeamProton;
    static const Type BeamNeutron;

    static const Particles_t& GetParticles() { return types; }

    class const_iterator : public Particles_t::const_iterator {
    public:

        const_iterator(const Particles_t::const_iterator& i):Particles_t::const_iterator(i)  {}

        const Type& operator*() const { return Particles_t::const_iterator::operator*().second; }

    };

    static const_iterator begin() { return const_iterator(types.begin()); }
    static const_iterator end()   { return const_iterator(types.end()); }

    static const TypeList_t& DetectableTypes() { return detectables; }
    static const TypeList_t& MCFinalStateTypes() { return mc_finalstate; }
    static const TypeList_t& NeutralMesons() { return neutral_mesons; }

    static const Type* GetTypeOfPlutoID(index_t pid);

    static mev_t CalculatePhotoproductionThreshold(mev_t m_sum, const Type& target);

};

inline bool operator<(const ant::ParticleTypeDatabase::Type& a, const ant::ParticleTypeDatabase::Type& b)
{
    return a.Name() < b.Name();
}

std::ostream& operator<< ( std::ostream& stream, const ant::ParticleTypeDatabase::Type& particle_type );

}



