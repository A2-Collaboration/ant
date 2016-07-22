#include "catch.hpp"

#include "analysis/input/goat/detail/TreeManager.h"
#include "analysis/input/goat/GoatReader.h"

#include "base/WrapTFile.h"

#include "TTree.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis::input;
void dotest_read();
void dotest_inputs();

TEST_CASE("GoatReader: Read ", "[analysis]") {
    dotest_read();
}

TEST_CASE("GoatReader: InputModules", "[analysis]") {
    dotest_inputs();
}


void dotest_read() {


    /// \todo Generate or read some Goat file?

//    auto filemanager = make_shared<WrapTFileInput>();
//    filemanager->OpenFile("NOTTHEREYET");
//    ant::analysis::input::GoatReader g(filemanager);

//    unsigned int n = 0;
//    while(g.hasData() && (n++ < 10)) {
//        auto event = g.ReadNextEvent();
//        cout << *event << endl;
//    }
}


class MyTreeRequestMgr: public TreeRequestManager {
protected:
    WrapTFileInput& m;
    TreeManager& t;

public:
    MyTreeRequestMgr(WrapTFileInput& _m, TreeManager& _t):
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

void dotest_inputs() {

    WrapTFileInput m;

    /// \todo Open some reasonable file here

    TreeManager treeManager;

    TrackInput tracks;
    TriggerInput trigger;
    DetectorHitInput detectorhit;
    ParticleInput photons("photons");
    TaggerInput tagger;

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

    bool init7 = tagger.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init7 == true) {
        cout << "Tagger Init OK" << endl;
    }


    if(init2 || init3 || init4 || init5 || init7) {

        cout << treeManager.GetEntries() << " entries " << endl;

        treeManager.GetEntry(1);
        if (init2) tracks.GetEntry();
        if (init3) trigger.GetEntry();
        if (init4) detectorhit.GetEntry();
        if (init5) photons.GetEntry();
        if (init7) tagger.GetEntry();

        if (init2) cout << "Goat Tracks: " << tracks.GetNTracks() << endl;
        if (init3) cout << "Trigger CBEsum: " << trigger.GetEnergySum() << "MeV" << endl;
        if (init4) cout << "DetectorHit: N NaIHits: " << detectorhit.NaI.Hits << endl;
        if (init5) cout << "Photons: " << photons.GetNParticles() << endl;
        if (init7) cout << "Tagger hits: " << tagger.GetNTagged() << endl;

    }

}