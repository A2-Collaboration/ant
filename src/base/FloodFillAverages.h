#pragma once

#include <vector>

namespace ant {

// suppose you have N elements
// getVal        = double(int i)
// setVal        = void(int i, double newval)
// getNeighbours = vector<int>(int i)
// getValid      = bool(int i)
template<typename GetVal_t, typename SetVal_t, typename GetNeighbours_t, typename GetValid_t>
inline void floodFillAverages(int N, GetVal_t getVal, SetVal_t setVal,
                              GetNeighbours_t getNeighbours, GetValid_t getValid) {


}

}