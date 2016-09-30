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
    static TGraph* getRatesVsTime(const std::list<treeLoader_t*>& tContainers, const analysis::HistogramFactory& histfac);
    static TGraph* getRatesVsTime(const std::list<treeLoader_t*>& tContainers, const size_t channel, const analysis::HistogramFactory& histfac);
    static TGraph* getLtVsTime(const std::list<treeLoader_t*>& tContainers, const analysis::HistogramFactory& histfac);
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
    static TID extractStartID(const std::string& f);

    const TID startID;

    analysis::HistogramFactory HistFac;

    treeLoader_t Bkg1;
    treeLoader_t Run;
    treeLoader_t Bkg2;



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

    taggEffTriple_t(const std::string& bkg1f, const std::string& runf, const std::string& bkg2f,
                    const analysis::HistogramFactory& histfac);

    std::string SetupName() const{return Bkg1.setupName;}

    void sanityChecks() const;

    const taggEff_t  GetTaggEffSubtracted() const;

    double GetDecayConstant() const;

};


}
}
}
