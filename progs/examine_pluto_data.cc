#include "analysis/input/goat/GoatReader.h"
#include "analysis/input/goat/detail/PlutoInput.h"
#include "base/Logger.h"

// Switch of some warnings for the Pluto headers
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#include "PParticle.h"
#pragma GCC diagnostic pop

#include <iostream>

using namespace std;
using namespace ant;

std::string plutoname(const int id) {
    if(id == 14001) return "gp";
    return makeStaticData()->GetParticleName(id);
}

int main(int argc, char** argv) {
    SetupLogger(argc, argv);

    input::GoatReader reader;

    for(int i=1; i<argc;++i) {
        if(argv[i][0] != '-')
            reader.AddInputFile(argv[i]);
    }

    reader.Initialize();

    const ant::input::PlutoInput& pluto = reader.GetPlutoInput();

    for(int n=0; n<100; ++n) {
    auto event = reader.ReadNextEvent();

    cout << "----------------------------" << endl;
    cout << pluto.Particles().size() << " Particles" << endl;

    int i=0;
    for(auto& p : pluto.Particles() ) {
        const PParticle* pl = p;
        cout << i++ << "\t" << plutoname(p->ID()) << "\t" << p->GetParentIndex() << "\t" << plutoname(pl->GetParentId()) << "\t|\t" <<  p->GetDaughterIndex() << "\t" <<  endl;
    }

    cout << "--- ant:" << endl;

    Particle::RecPrint(event->MCTrue().Intermediates().Get(ParticleTypeDatabase::BeamProton).front());
    cout << endl;

    if(!reader.hasData())
        break;
    }


    return 0;

}
