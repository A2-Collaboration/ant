#pragma once
#include <string>

class TTree;

namespace ant {

struct KinfitExtract {
    static void Sweep1(TTree* t);

    static std::string JSONTest();
    static std::string LoadAndDump(const std::string& s);

    static void short_string_test();
};


}
