#ifndef EXTRACT_RESOLUTIONS_H
#define EXTRACT_RESOLUTIONS_H

#include <string>

class TTree;
class TCut;
class TF1;



class ExtractResolutions {
public:

    static void AnalyseThetaCB(TTree* tree);
    static void AnalysePhiCB(TTree* tree);
    static void AnalyseECB(TTree* tree);

    static TF1* LogNormal();
    static TF1* SigmaFit();
};

#endif
