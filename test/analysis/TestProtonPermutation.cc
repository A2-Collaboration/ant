/**
  * @brief Test for ProtonPermutation.
  *   Use Candidate cluster sizes to track candidates and particles
  */

#include "catch.hpp"
#include "catch_config.h"

#include "analysis/utils/ProtonPermutation.h"

#include "tree/TCandidate.h"
#include <memory>
#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::analysis;


/**
 * @brief Create a candidate list, cluster sizes as given in vector
 * @param z
 * @return
 */
TCandidatePtrList MakeCands(const std::vector<u_int16_t>& z) {
    TCandidatePtrList l;

    for(const auto& i : z) {
        auto c = std::make_shared<TCandidate>();
        c->ClusterSize = i;
        l.emplace_back(move(c));
    }

    return l;
}

void dotest_permutation(const TCandidatePtrList& cands, const vector<u_int16_t> &numbers);

TEST_CASE("Analysis: Proton Permutation", "[analysis]") {
    const vector<u_int16_t> numbers = {0,1,2,3,4};
    const auto cands = MakeCands(numbers);
    dotest_permutation(cands, numbers);
}

template <typename T>
void sort(std::vector<T>& v) {
    std::sort(v.begin(), v.end());
}

void dotest_permutation(const TCandidatePtrList& cands, const vector<u_int16_t>& numbers) {

    assert(numbers.size() >= 1);



    unsigned n=0;

    vector<u_int16_t> a;
    for(utils::ProtonPermutation perm(cands); perm.Good(); perm.Next()) {

        vector<u_int16_t> b;

        REQUIRE(perm.Proton());
        REQUIRE(perm.Proton()->Type() == ParticleTypeDatabase::Proton);
        a.push_back(perm.Proton()->Candidate->ClusterSize);
        b.push_back(perm.Proton()->Candidate->ClusterSize);

        REQUIRE(perm.Photons().size() == numbers.size()-1);
        for(const auto& p : perm.Photons()) {
            REQUIRE(p->Type() == ParticleTypeDatabase::Photon);
            b.push_back(p->Candidate->ClusterSize);
        }

        sort(b);

        REQUIRE(b == numbers); // make sure all candidates are accounted for (either as proton or in the photons list)

        ++n;

    }

    sort(a);

    REQUIRE(a == numbers);  // make sure every candidate has been the proton

    REQUIRE(n==cands.size());
}
