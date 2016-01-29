#ifndef OMEGAETAG_H
#define OMEGAETAG_H

class TTree;
class TCut;
class TFile;

class OmegaEtaG {
public:

    static void Chi2Plot(TTree* tree);

    static void Chi2CutCompare(TTree* tree, double chi2=4);

    static void gggIM_MM(TTree* tree, const TCut& cut);
    static void ggIM(TTree* tree, const TCut& cut);

    static void Analyse(TTree* tree);

    static void kinfit1();

    static void Plot(TTree* tree, const double binscale=1.0);
    static void DataMC(TFile* mc_file, TFile* data_file, const double mcscale=1.0);
};

#endif
