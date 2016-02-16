#ifndef OMEGAETAG_H
#define OMEGAETAG_H

class TTree;
class TCut;
class TFile;

class OmegaEtaG {
public:

    static void Plot(TTree* tree, const TCut& extra_cut, const double binscale=1.0);
    static void PlotBGs(TTree* tree, const double binscale=1.0);
    static void DataMC(TFile* mc_file, TFile* data_file, const double mcscale=1.0);
    static void DataMCBGs(TFile* mc_file, TFile* data_file, const double mcscale);
};

#endif
