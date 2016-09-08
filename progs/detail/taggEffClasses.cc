
#include "analysis/physics/common/ProcessTaggEff.h"
#include "analysis/plot/HistogramFactories.h"

#include "base/Logger.h"
#include "base/WrapTFile.h"

#include "expconfig/detectors/EPT.h"

#include "tree/TAntHeader.h"



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
    taggEff_t(const string& setup, const TID& firstID):
        Setup(setup),
        FirstID(firstID),
        TaggEffs(),
        TaggEffErrors() {}
};

struct treeContainer_t
{
    WrapTFileInput wrapFile;
    physics::ProcessTaggEff::TreeScalarReads wrapTree;

    string setupName;
    size_t nchannels;


    treeContainer_t(const string& filename):
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

    vector<pair<double,double>> getLiveTimes()
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

class taggEffTriple_t
{
protected:

    treeContainer_t bkg1;
    treeContainer_t run;
    treeContainer_t bkg2;

    TID startID;


public:

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
    }

    string SetupName() const{return bkg1.setupName;}

   unsigned sanityChecks(const treeContainer_t::means_t& bkg1,
                 const treeContainer_t::means_t& run,
                 const treeContainer_t::means_t& bkg2) const
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



    const taggEff_t  GetTaggEff() const
    {
        taggEff_t result(SetupName(),startID);

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
