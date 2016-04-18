#pragma once
#include <string>

class TTree;
class TH1;
class TDirectory;
class THStack;

namespace ant {

struct KinfitExtract {
    static void Sweep1(TTree* t);

    static std::string JSONTest();
    static std::string LoadAndDump(const std::string& s);

    static void short_string_test();


    static THStack* DrawPullSame(TH1* data, TH1* mc, const std::string& title="", const bool draw=true);
    static void DrawPullSameAll(TDirectory* data, TDirectory* mc);
};


}
