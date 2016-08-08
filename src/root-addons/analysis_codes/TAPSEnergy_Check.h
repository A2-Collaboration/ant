#pragma once

class TFile;

namespace ant {

struct TAPSEnergy_Check {

    static void Analyse(TFile* file);

    static void AnalyseTree(TFile* file);

};

}