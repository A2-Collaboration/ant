#include <vector>
#include <memory>
#include <iostream>
#include <string>
#include <forward_list>

#include "tree/MemoryPool.h"
#include "tree/TDetectorRead.h"

using namespace std;
using namespace ant;



int main()
{
    vector<MemoryPool<TDetectorRead>::Item> items;
    for(int i=0;i<10;i++) {
        items.push_back(MemoryPool<TDetectorRead>::Get());
    }
    items.clear();

    auto item = MemoryPool<TDetectorRead>::Get();
//    item = MemoryPool<TDetectorRead>::Get();
//    item = MemoryPool<TDetectorRead>::Get();
//    item = MemoryPool<TDetectorRead>::Get();
//    item = MemoryPool<TDetectorRead>::Get();
//    item = MemoryPool<TDetectorRead>::Get();


    return 0;
}
