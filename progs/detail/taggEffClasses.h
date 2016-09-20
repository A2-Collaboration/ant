#pragma once
#include <string>
#include <vector>
#include <list>


#include "tree/TID.h"
#include "base/WrapTFile.h"

#include "analysis/physics/common/ProcessTaggEff.h"


class TGraph;

namespace ant
{
namespace progs
{
namespace taggeff
{

struct treeLoader_t;
struct timedData
{
    static TGraph* getRatesVsTime(const std::list<treeLoader_t*>& tContainers);
    static TGraph* getRatesVsTime(const std::list<treeLoader_t*>& tContainers, const size_t channel);
    static TGraph* getLtVsTime(const std::list<treeLoader_t*>& tContainers);
};



struct taggEff_t
{
    std::string          Setup;
    TID             FirstID;
    std::vector<double>  TaggEffs;
    std::vector<double>  TaggEffErrors;
    std::vector<double>  BkgFitChi2;
    taggEff_t(const std::string& setup, const TID& firstID, const size_t nChannels):
        Setup(setup),
        FirstID(firstID),
        TaggEffs(nChannels),
        TaggEffErrors(nChannels),
        BkgFitChi2(nChannels){}
};


struct treeLoader_t
{
    WrapTFileInput wrapFile;
    ant::analysis::physics::ProcessTaggEff::TreeScalarReads wrapTree;

    std::string setupName;
    size_t nchannels;

    uint32_t startTime;


    treeLoader_t(const std::string& filename);

    struct means_t {
        std::vector<double> Scalers;
        std::vector<double> ScalersErrors;
        double Livetime = 0;
        double LivetimeError = 0;
        std::vector<double> Tdcs;
        std::vector<double> TdcsErrors;
    };

    means_t getMeans() const;

    TTree* Tree() const { return wrapTree.Tree; }

    std::vector<std::pair<double,double>> getLiveTimes() const;
};

class taggEffTriple_t
{
protected:

    treeLoader_t Bkg1;
    treeLoader_t Run;
    treeLoader_t Bkg2;


    TID startID;

    void initBkgFits();

public:

    TGraph* AvgBkgRates;
    TF1*    AvgBkgFit;
    struct bkgFit_t
    {
        TGraph* Graph;
        TF1*    Fit;

        bkgFit_t(TGraph* graph, const size_t channel);

        void doFit(const IntervalD& fitrange, const double lambda);
        double operator ()(const double time)  const;
    };

    std::vector<bkgFit_t> bkgFits;



    TGraph* avgRatesSub = nullptr;
    TGraph* avgRates    = nullptr;

    taggEffTriple_t(const std::string& bkg1f, const std::string& runf, const std::string& bkg2f);

    std::string SetupName() const{return Bkg1.setupName;}

    unsigned sanityChecks(const treeLoader_t::means_t& bkg1,
                          const treeLoader_t::means_t& run,
                          const treeLoader_t::means_t& bkg2) const;

    const taggEff_t  GetTaggEffSubtracted() const;

    const taggEff_t  GetTaggEff() const;

};


}
}
}
