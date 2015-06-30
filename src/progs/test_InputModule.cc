#include "input/FileManager.h"
#include "input/TreeManager.h"
#include "input/PlutoInput.h"
#include "input/TrackInput.h"
#include "input/TriggerInput.h"
#include "input/DetectorHitInput.h"
#include "input/ParticleInput.h"
#include "input/GeantInput.h"
#include "input/TaggerInput.h"

#include "TTree.h"

#include <iostream>

using namespace ant::input;
using namespace std;

class MyTreeRequestMgr: public TreeRequestManager {
protected:
    FileManager& m;
    TreeManager& t;

public:
    MyTreeRequestMgr(FileManager& _m, TreeManager& _t):
        m(_m), t(_t) {}

    TTree *GetTree(const string &name) {
        TTree* tree = nullptr;
        if( m.GetObject(name, tree) ) {
            t.AddTree(tree);
            return tree;
        } else
            return nullptr;
    }
};

int main(int argc, char** argv) {

    FileManager m;

    for(int i = 1; i < argc; ++i) {
        m.OpenFile(argv[i]);
    }

    TreeManager treeManager;

    PlutoInput pluto;
    TrackInput tracks;
    TriggerInput trigger;
    DetectorHitInput detectorhit;
    ParticleInput photons("photons");
    GeantInput geant;
    TaggerInput tagger;

    bool init1 = pluto.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init1 == true) {
        cout << "Pluto Init OK" << endl;
    }

    bool init2 = tracks.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init2 == true) {
        cout << "Tracks Init OK" << endl;
    }

    bool init3 = trigger.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init3 == true) {
        cout << "Trigger Init OK" << endl;
    }

    bool init4 = detectorhit.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init4 == true) {
        cout << "DetectorHit Init OK" << endl;
    }

    bool init5 = photons.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init5 == true) {
        cout << "Photons Init OK" << endl;
    }

    bool init6 = geant.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init6 == true) {
        cout << "Geant Init OK" << endl;
    }

    bool init7 = tagger.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init7 == true) {
        cout << "Tagger Init OK" << endl;
    }


    if( init1 || init2 || init3 || init4 || init5 || init6 || init7) {

        cout << treeManager.GetEntries() << " entries " << endl;

        treeManager.GetEntry(1);
        if (init1) pluto.GetEntry();
        if (init2) tracks.GetEntry();
        if (init3) trigger.GetEntry();
        if (init4) detectorhit.GetEntry();
        if (init5) photons.GetEntry();
        if (init6) geant.GetEntry();
        if (init7) tagger.GetEntry();

        if (init1) cout << "Pluto Particles: " << pluto.Particles().size() << endl;
        if (init2) cout << "Goat Tracks: " << tracks.GetNTracks() << endl;
        if (init3) cout << "Trigger CBEsum: " << trigger.GetEnergySum() << "MeV" << endl;
        if (init4) cout << "DetectorHit: N NaIHits: " << detectorhit.GetNNaIHits() << endl;
        if (init5) cout << "Photons: " << photons.GetNParticles() << endl;
        if (init6) cout << "Geant: " << geant.GetNCBHits() << endl;
        if (init7) cout << "Tagger hits: " << tagger.GetNTagged() << endl;

    }


    m.CloseAll();

}
