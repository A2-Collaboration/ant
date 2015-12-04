#ifndef SYMMETRIC_PI0_H
#define SYMMETRIC_PI0_H

class TTree;
class TH2;

class SymmetricPi0 {
public:

    static void Analyse(TTree* tree);
    static void findMax(TH2* h);
};


#endif
