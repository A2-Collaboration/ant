#ifndef OMEGAETAG_H
#define OMEGAETAG_H

class TTree;
class TCut;

class OmegaEtaG {
public:

    static void Chi2Plot(TTree* tree);

    static void Chi2CutCompare(TTree* tree, double chi2=4);

    static void gggIM_MM(TTree* tree, const TCut& cut);

    static void Analyse(TTree* tree);

    static void kinfit1();
};

#endif
