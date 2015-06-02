#include <iostream>
#include <typeinfo>
#include <vector>

#include "base/types.h"
#include "ParticleType.h"
#include "Track.h"
#include "Particle.h"

#include "utils/combinatorics.h"

#include "base/interval.h"

#include "Detector.h"
#include "utils/matcher.h"

using namespace std;
using namespace ant;

int main() {

    ParticleTypeDatabase::Print();

    for( ParticleTypeDatabase::const_iterator a = ParticleTypeDatabase::begin(); a!=ParticleTypeDatabase::end(); a++) {
        cout << *a << endl;
    }

    for( auto& a : ParticleTypeDatabase() ) {
        cout << endl;
        for( auto& b : ParticleTypeDatabase() ) {
            cout << a << " " << b << "    " << (a == b) << endl;
        }
    }

    TrackPtr t(new Track(mev_t(100),
            radian_t(2.0),
            radian_t(1.0),
            time_t(0.324),
            clustersize_t(4),
            detector_t::NaI,
            mev_t(0.4),
            mev_t(0.3),
            mev_t(0.2)
            ));
    cout << t << endl;

    Particle g(ParticleTypeDatabase::Photon, t);
    cout << g << endl;
    g.ChangeType(ParticleTypeDatabase::Proton);
    cout << g << endl;
    g.ChangeType(ParticleTypeDatabase::PiCharged);
    cout << g << endl;

    cout << "=== combinatorics test ===" <<endl;

    vector<int> numbers = {1,2,3,4,5};

    KofNvector<int> combinations(numbers,4);
    do {

        for( auto& i: combinations ) {
            cout << i << " ";
        }

        cout << endl;

    } while(combinations.next());

    cout << " as for loop " << endl;
    for(auto ddd = makeCombination(numbers,4); !ddd.Done(); ++ddd ) {

        for( auto& i: ddd ) {
            cout << i << " ";
        }
        cout << endl;
    }

    cout << "----" << endl;
    vector<int> numbers2 = {1};

    KofNvector<int> combinations2(numbers2,2);
    do {

        for( auto& i: combinations2 ) {
            cout << i << " ";
        }

        cout << endl;

    } while(combinations2.next());

    interval<double> i;
    cout << i << endl;

    i.SetWidth(5);
    cout << i << endl;

    i.SetCenter(2.5);
    cout << i << endl;

    cout << i.Center() << " " << i.Length() << endl;

    cout << "========================" << endl;

    detector_t d = detector_t::NaI;

    if( d & detector_t::anyCB )
        cout << "is in cb";
    else
        cout << "NOT" << endl;

    cout << d << endl;

    d ^= detector_t::PbWO4;
    cout << d << endl;

    std::string detector_name = detector_t::anyTAPS;

    cout << detector_name << endl;

    cout << "---------------------" <<endl;

    ParticleList mctrue;
    mctrue.emplace_back( new Particle( ParticleTypeDatabase::Photon, 100, 1,2));
    mctrue.emplace_back( new Particle( ParticleTypeDatabase::Photon, 200, .2 ,.4));

    ParticleList rec;

    TrackPtr r1(new Track(mev_t(100),
            radian_t(1.01),
            radian_t(2.01),
            time_t(0.324),
            clustersize_t(4),
            detector_t::NaI,
            mev_t(0.4),
            mev_t(0.3),
            mev_t(0.2)
            ));
    rec.emplace_back( new Particle(ParticleTypeDatabase::Photon, r1));

    TrackPtr r2(new Track(mev_t(210),
            radian_t(.19),
            radian_t(.41),
            time_t(0.323),
            clustersize_t(5),
            detector_t::NaI,
            mev_t(0.4),
            mev_t(0.3),
            mev_t(0.2)
            ));

    rec.emplace_back( new Particle(ParticleTypeDatabase::Photon, r2));

    TrackPtr r3(new Track(mev_t(210),
            radian_t(1.19),
            radian_t(.41),
            time_t(0.323),
            clustersize_t(5),
            detector_t::NaI,
            mev_t(0.4),
            mev_t(0.3),
            mev_t(0.2)
            ));

    rec.emplace_back( new Particle(ParticleTypeDatabase::Photon, r3));

    // match by angle
    //auto matched = Match(mctrue,rec, matchAngle);

    // match by energy using labda function
    auto matched = utils::match1to1(mctrue,rec, [] (const ParticlePtr& a, const ParticlePtr& b) {return abs(a->E() - b->E());});

    cout << "matches: " << matched.size() << endl;

    for( auto& match : matched) {
        cout << "\n" << match.a << "\n" << match.b << endl;
        cout << "E diff="<< abs(match.a->E() - match.b->E()) << " score=" << match.score << endl;
    }

    cout << " or with numbers: Only accept if distance < 3:" << endl;

    std::list<int> n1 = { 100,1, 10 , 20 };
    std::vector<int> n2 = { 15, 30 , 8, 5, 99 };

    auto imatched = utils::match1to1(n1, n2, utils::matchDistance<int>, IntervalD(0,3));

    for( auto& match : imatched) {
        cout << "\n" << match.a << " <-> " << match.b << "  score: " << match.score << endl;
    }

    cout << "===========================" << endl;

    return 0;

}
