#pragma once

#include <vector>
#include <limits>
#include <algorithm>

#include "base/std_ext/container.h"

namespace ant {

// suppose you have N elements
// getVal        = double(int i)
// setVal        = void(int i, double newval)
// getNeighbours = vector<int>(int i)
// getValid      = bool(int i)
template<typename GetVal_t, typename SetVal_t, typename GetNeighbours_t, typename GetValid_t>
inline void floodFillAverages(int N, GetVal_t getVal, SetVal_t setVal,
                              GetNeighbours_t getNeighbours, GetValid_t getValid)
{
    struct invalid_t {
        int Index;
        explicit invalid_t(int i) : Index(i) {}
        int ValidNeighbours = 0;
        std::vector<int> InvalidNeighbours;
        double Value = std::numeric_limits<double>::quiet_NaN();
        bool Visited = false;
    };

    std::vector<invalid_t> invalids;

    for(int i=0;i<N;i++) {
        if(getValid(i))
            continue;
        invalids.emplace_back(i);
        for(auto j : getNeighbours(i)) {
            if(getValid(j))
                invalids.back().ValidNeighbours++;
            else
                invalids.back().InvalidNeighbours.push_back(j);
        }
    }

    // average calculation takes into account already set invalids
    auto getAvg = [&invalids, getNeighbours, getVal, getValid] (int i) {
        double sum = 0;
        int n = 0;
        // average over valid neighbours
        for(int j : getNeighbours(i)) {
            if(getValid(j)) {
                sum += getVal(j);
                n++;
            }
            else {
                // and possibly visited invalids
                auto it = std::find_if(invalids.begin(), invalids.end(), [j] (const invalid_t& i) {
                    return j == i.Index && i.Visited;
                });
                if(it != invalids.end()) {
                    sum += it->Value;
                    n++;
                }
            }
        }
        return sum/n;
    };

    while(true) {
        std::sort(invalids.begin(), invalids.end(),
                  [] (const invalid_t& a, const invalid_t& b) {
            return a.ValidNeighbours > b.ValidNeighbours;
        });
        // find the first unvisited item (if any left)
        // with the highest number of ValidNeighbours
        auto it = std::find_if(invalids.begin(), invalids.end(), [] (const invalid_t& i) {
            return !i.Visited;
        });
        if(it == invalids.end())
            break;

        // there might be more than one unvisited invalid
        // with same number of valid neighbours
        std::vector<decltype(it)> unvisited{it++};
        for(; it != invalids.end(); ++it) {
            if(it->ValidNeighbours < unvisited.front()->ValidNeighbours)
                break; // stop immediately as invalids are sorted
            if(it->Visited)
                continue; // skip visited
            unvisited.emplace_back(it);
        }

        for(const auto& it : unvisited) {

            // important to set the value for getAvg in next iteration
            it->Value = getAvg(it->Index);
            setVal(it->Index, it->Value);

            // neighbour relatings might not be reflexive
            // so we cannot iterate over getNeighbours(i.Index),
            // but need to check if it->Index is neighbour of all other invalids
            for(auto& invalid : invalids) {
                if(invalid.Index == it->Index)
                    continue;

                if(std_ext::contains(invalid.InvalidNeighbours, it->Index)) {
                    invalid.ValidNeighbours++;
                }
            }

            it->Visited = true;

        }
    }
}

}