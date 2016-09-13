
#include "analysis/physics/common/ProcessTaggEff.h"
#include "analysis/plot/HistogramFactories.h"

#include "base/Logger.h"
#include "base/WrapTFile.h"

#include "expconfig/detectors/EPT.h"

#include "tree/TAntHeader.h"
#include "TF1.h"
#include "TGraph.h"


using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;



struct taggEff_t
{
    string          Setup;
    TID             FirstID;
    vector<double>  TaggEffs;
    vector<double>  TaggEffErrors;
    taggEff_t(const string& setup, const TID& firstID, const size_t nChannels):
        Setup(setup),
        FirstID(firstID),
        TaggEffs(nChannels),
        TaggEffErrors(nChannels) {}
};

static size_t FillGraph(TGraph* graph, const double x, const double y)
{
    graph->SetPoint(graph->GetN(),x,y);
    return graph->GetN();
}

struct treeLoader_t
{
    WrapTFileInput wrapFile;
    physics::ProcessTaggEff::TreeScalarReads wrapTree;

    string setupName;
    size_t nchannels;

    uint32_t startTime;


    treeLoader_t(const string& filename):
        wrapFile(filename)
    {
        auto treeName = physics::ProcessTaggEff::treeAccessName();
        if (!wrapFile.GetObject(treeName,wrapTree.Tree))
        {
            LOG(ERROR) << "Cannot find tree " << treeName
                       << " in file " << filename;
            exit(EXIT_FAILURE);
        }

        wrapTree.LinkBranches(wrapTree.Tree);

        ant::TAntHeader* h;
        wrapFile.GetObject("AntHeader",h);
        setupName = h->SetupName;
        ExpConfig::Setup::SetManualName(setupName);
        nchannels = ExpConfig::Setup::GetDetector<TaggerDetector_t>()->GetNChannels();

        wrapTree.Tree->GetEntry(0);
        startTime = wrapTree.EvID.Value->Timestamp;
    }

    struct means_t {
        vector<double> Scalers;
        vector<double> ScalersErrors;
        double Livetime = 0;
        double LivetimeError = 0;
        vector<double> Tdcs;
        vector<double> TdcsErrors;
    };

    means_t getMeans() const
    {
        means_t means;
        means.Scalers.resize(nchannels);
        means.ScalersErrors.resize(nchannels);
        means.Tdcs.resize(nchannels);
        means.TdcsErrors.resize(nchannels);

        vector<std_ext::RMS> rms_scalers(nchannels);
        vector<std_ext::RMS> rms_tds(nchannels);
        std_ext::RMS rms_livetime;


        auto nentries = wrapTree.Tree->GetEntries();

        for (auto entry = 0 ; entry < nentries ; ++entry)
        {
            wrapTree.Tree->GetEntry(entry);
            for ( auto channel = 0u ; channel < nchannels; ++channel)
            {
                rms_scalers.at(channel).Add(1.0 * wrapTree.TaggRates().at(channel));
                rms_tds.at(channel).Add(1.0 * wrapTree.TDCRates().at(channel));
            }
            rms_livetime.Add(1.0 * wrapTree.ExpLivetime);
        }

        for ( auto channel = 0u ; channel < nchannels ; ++channel)
        {
            means.Scalers.at(channel)       = rms_scalers.at(channel).GetMean();
            means.ScalersErrors.at(channel) = rms_scalers.at(channel).GetSigmaMean();

            means.Tdcs.at(channel)          = rms_tds.at(channel).GetMean();
            means.TdcsErrors.at(channel)    = rms_tds.at(channel).GetSigmaMean();
        }
        means.Livetime      = rms_livetime.GetMean();
        means.LivetimeError = rms_livetime.GetSigmaMean();

        return means;
    }

    TTree* Tree() const { return wrapTree.Tree; }

    vector<pair<double,double>> getLiveTimes() const
    {
        vector<pair<double,double>> vTimeLivetime;

        auto nentries = wrapTree.Tree->GetEntries();
        for ( auto en = 0u ; en < nentries ; ++en )
        {
            wrapTree.Tree->GetEntry(en);
            vTimeLivetime.emplace_back(pair<double,double>({wrapTree.Clock,wrapTree.ExpLivetime}));
        }

        return vTimeLivetime;
    }


};

static TGraph* getRatesVsTime(const list<treeLoader_t*>& tContainers)
{
    auto graph2d = new TGraph();

    //in vain for our canvas, but ...
    graph2d->SetMarkerStyle(kPlus);
    graph2d->GetXaxis()->SetTitle("time [s]");
    graph2d->GetYaxis()->SetTitle("avg. rate [Hz]");

    auto first_time = numeric_limits<uint32_t>::quiet_NaN();
    double timeInRun(0);

    for ( auto t : tContainers )
    {
        timeInRun = 0;
        for ( auto en = 0u ; en < t->Tree()->GetEntries() ; ++en)
        {
            t->Tree()->GetEntry(en);
            if (first_time == numeric_limits<uint32_t>::quiet_NaN() && en == 0)
                first_time = t->wrapTree.EvID.Value->Timestamp;


            timeInRun += t->wrapTree.Clock() / 1.0e6;

            auto evTime = (  timeInRun
                             + t->wrapTree.EvID.Value->Timestamp
                             - first_time );

            std_ext::RMS rmsRate;
            for ( const auto& tr: t->wrapTree.TaggRates())
                rmsRate.Add(tr);

            FillGraph(graph2d,evTime, rmsRate.GetMean());
        }
    }

    return graph2d;
}

class taggEffTriple_t
{
protected:

    treeLoader_t bkg1;
    treeLoader_t run;
    treeLoader_t bkg2;


    TID startID;

    void initBkgFuntion()
    {
        bkgGraph = getRatesVsTime({addressof(bkg1),addressof(bkg2)});

        double tmax(0);
        double dummy(0);
        bkgGraph->GetPoint(bkgGraph->GetN()-1,tmax,dummy);
        bkgFit = new TF1("bkgFit","[0] + [1] * exp( - [2] * x)",
                         0, tmax);
        bkgGraph->Fit("bkgFit");

        LOG(INFO) << " bkg-fit: b(t) = "
                  << bkgFit->GetParameter(0) << " + "
                  << bkgFit->GetParameter(1) << " * exp( - "
                  << bkgFit->GetParameter(2) << " * t )";
    }

public:

    TF1*    bkgFit      = nullptr;
    TGraph* bkgGraph    = nullptr;
    TGraph* avgRates    = nullptr;
    TGraph* avgRatesSub = nullptr;

    taggEffTriple_t(const string& bkg1f_, const string& runf_, const string& bkg2f_):
        bkg1(bkg1f_),
        run(runf_),
        bkg2(bkg2f_)
    {
        LOG(INFO) << "Loading TaggEff measurement for:  "
                  << "1st bkg: " << bkg1f_ << "  "
                  << "run: "     << runf_  << "  "
                  << "2nd bkg: " << bkg2f_;

        if ( !(bkg1.setupName == bkg2.setupName &&
             bkg2.setupName == run.setupName ) )
            throw runtime_error("Files in TaggEff-triple not from same Setup!");
        WrapTFileInput file(bkg2f_);
        ant::TAntHeader* header;
        file.GetObject("AntHeader",header);
        if (!header)
            throw runtime_error("No Ant header in found!");
        startID = header->LastID;

        initBkgFuntion();
        avgRates = new TGraph();
        avgRates->SetMarkerStyle(kPlus);
        avgRates->SetMarkerColor(kBlue);
        avgRatesSub = new TGraph();
        avgRatesSub->SetMarkerStyle(kPlus);
        avgRatesSub->SetMarkerColor(kRed);
    }

    string SetupName() const{return bkg1.setupName;}

   unsigned sanityChecks(const treeLoader_t::means_t& bkg1,
                 const treeLoader_t::means_t& run,
                 const treeLoader_t::means_t& bkg2) const
    {
        auto severity = 0u;
        auto nchannels = bkg1.Scalers.size();




        // per file checks:
        for ( auto m: {bkg1, run, bkg2})
        {
            for (auto ch = 0u ; ch < nchannels ; ++ch)
            {
                if (m.Scalers.at(ch) < m.Tdcs.at(ch))
                {
                    LOG(ERROR) << "Detected file with higher TDC rate than scaler rate!";
                    LOG(ERROR) << m.Scalers.at(ch) << " < " << m.Tdcs.at(ch);
                    severity++;
                }
            }
            if ( !(0 <= m.Livetime && m.Livetime <= 1))
            {
                LOG(ERROR) << "Detected file with live time not in [0,1]!";
                LOG(ERROR) << "livetime: " << m.Livetime;
                severity++;
            }
        }

        // bkg
        for ( auto ch = 0u ; ch < nchannels ; ++ch)
        {
            if ( run.Scalers.at(ch) < bkg1.Scalers.at(ch) || run.Scalers.at(ch) < bkg2.Scalers.at(ch))
                LOG(WARNING) << "Detected file with higher background than actual measurement in scalers, channel " << ch << "!";
            if ( run.Tdcs.at(ch) < bkg1.Tdcs.at(ch) || run.Tdcs.at(ch) < bkg2.Tdcs.at(ch))
                LOG(WARNING) << "Detected file with higher background than actual measurement in TDCs, channel " << ch << "!";
        }
        return severity;
    }

   const taggEff_t  GetTaggEffSubtracted() const
   {
       taggEff_t result(SetupName(),startID,run.nchannels);

       auto formula = [] (const double TDC,
                          const double L,
                          const double scaler)
       {
           return 1.0 * TDC / (L * scaler);
       };

       auto formula_err = [] (const double err_a, const double TDC,
                          const double a,
                          const double b)
       {
           return 1.0 * err_a * TDC / (a*a*b);
       };

       vector<std_ext::RMS> TDCs(run.nchannels);
       vector<std_ext::RMS> scalers(run.nchannels);
       std_ext::RMS         L;

       auto time = 1.0 * run.startTime - bkg1.startTime;
       auto nentries = run.wrapTree.Tree->GetEntries();

       for (auto entry = 0 ; entry < nentries ; ++entry)
       {
           run.wrapTree.Tree->GetEntry(entry);
           time += run.wrapTree.Clock() / 1.0e6;
           std_ext::RMS graphTDC;
           std_ext::RMS graphTDCsub;
           for (auto  ch = 0u ; ch < run.nchannels ; ++ch)
           {
               TDCs.at(ch).Add(run.wrapTree.TDCRates().at(ch));
               auto rate_sub = run.wrapTree.TaggRates().at(ch) - bkgFit->Eval(time);
               scalers.at(ch).Add( rate_sub );
               graphTDC.Add(run.wrapTree.TaggRates().at(ch));
               graphTDCsub.Add(rate_sub);
           }
           FillGraph(avgRates,time,graphTDC.GetMean());
           FillGraph(avgRatesSub,time,graphTDCsub.GetMean());
           L.Add(run.wrapTree.ExpLivetime);
       }

       for ( auto ch = 0u ; ch < run.nchannels ; ++ch)
       {
           result.TaggEffs.at(ch) = formula(TDCs.at(ch).GetMean(),
                                            L.GetMean(),
                                            scalers.at(ch).GetMean() );
           result.TaggEffErrors.at(ch) = sqrt(  sqr(    formula(TDCs.at(ch).GetSigmaMean(), L.GetMean(), scalers.at(ch).GetMean()))
                                              + sqr(formula_err(L.GetSigmaMean(),TDCs.at(ch).GetMean(),
                                                                L.GetMean(),scalers.at(ch).GetMean()))
                                              + sqr(formula_err(scalers.at(ch).GetSigmaMean(),TDCs.at(ch).GetMean(),
                                                                scalers.at(ch).GetMean(),L.GetMean()))
                                             );
       }


       return result;
   }

    const taggEff_t  GetTaggEff() const
    {
        taggEff_t result(SetupName(),startID,run.nchannels);

        auto denum = [] (const double L,  const double s,
                          const double L1, const double s1,
                          const double L2, const double s2)
        {
            return L * s - ( ( L1 * s1 + L2 * s2 ) / 2.0 );
        };
        auto denum_bkg = [] (const double L,  const double s,
                          const double L1, const double s1,
                          const double L2, const double s2)
        {
            return 2 * L *s - L1 * s1 - L2 * s2;
        };


        const auto m_bkg1 = bkg1.getMeans();
        const auto m_run = run.getMeans();
        const auto m_bkg2 = bkg2.getMeans();
        if ( sanityChecks(m_bkg1,m_run,m_bkg2) > 0 )
            throw runtime_error("Tagg-Eff-triple didn't pass sanity checks, doublecheck input files!");

        auto nchannels = bkg1.nchannels;
        result.TaggEffs.resize(nchannels,0);
        result.TaggEffErrors.resize(nchannels,0);

        for (auto ch = 0u ; ch < nchannels ; ++ch)
        {
            auto denum_ch         = denum(m_run.Livetime,  m_run.Scalers.at(ch),
                                          m_bkg1.Livetime, m_bkg1.Scalers.at(ch),
                                          m_bkg2.Livetime, m_bkg2.Scalers.at(ch));
            auto denum_bkg_ch_sqr = denum_bkg(m_run.Livetime,  m_run.Scalers.at(ch),
                                          m_bkg1.Livetime, m_bkg1.Scalers.at(ch),
                                          m_bkg2.Livetime, m_bkg2.Scalers.at(ch));


            result.TaggEffs.at(ch)  = 1.0 * m_run.Tdcs.at(ch) / denum_ch;

            result.TaggEffErrors.at(ch) =
                    sqrt(   sqr( m_run.TdcsErrors.at(ch)/denum_ch )
                          + sqr( m_run.LivetimeError * m_run.Tdcs.at(ch) * m_run.Scalers.at(ch)  / sqr(denum_ch) )
                          + sqr( m_run.ScalersErrors.at(ch) * m_run.Livetime * m_run.Tdcs.at(ch) / sqr(denum_ch) )
                          + sqr( m_bkg1.LivetimeError * 2.0 * m_run.Tdcs.at(ch) * m_bkg1.Scalers.at(ch)  / denum_bkg_ch_sqr )
                          + sqr( m_bkg1.ScalersErrors.at(ch) * 2.0 * m_run.Tdcs.at(ch) * m_bkg1.Livetime / denum_bkg_ch_sqr )
                          + sqr( m_bkg2.LivetimeError * 2.0 * m_run.Tdcs.at(ch) * m_bkg2.Scalers.at(ch)  / denum_bkg_ch_sqr )
                          + sqr( m_bkg2.ScalersErrors.at(ch) * 2.0 * m_run.Tdcs.at(ch) * m_bkg2.Livetime / denum_bkg_ch_sqr )
                         );

            if ( !(IntervalD(0,1).Contains(result.TaggEffs.at(ch))) )
                LOG(WARNING)  << "Calculated Tagging efficiency for channel " << ch
                              << " not in [0,1]: "
                              << result.TaggEffs.at(ch) << "!";
        }

        return result;
    }

};
