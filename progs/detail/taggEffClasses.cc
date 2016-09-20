#include "taggEffClasses.h"

#include "analysis/plot/HistogramFactories.h"


#include "base/Logger.h"
#include "base/interval.h"

#include "expconfig/detectors/EPT.h"

#include "tree/TAntHeader.h"
#include "TF1.h"
#include "TGraph.h"

#include "tools.h"

using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;
using namespace ant::progs::taggeff;
using namespace ant::progs::tools;

const vector<double> startparams({0.66,1.60,0.001});
constexpr double upperlimit = 100;

template<typename Func>
void runOverContainer(const list<treeLoader_t*>& tContainers, Func f) {

    uint32_t first_time = 0;
    bool first_time_valid = false;

    for(const auto& t : tContainers )
    {
        double timeInRun = 0;
        for ( auto en = 0u ; en < t->Tree()->GetEntries() ; ++en)
        {
            t->Tree()->GetEntry(en);
            if (!first_time_valid && en == 0) {
                first_time = t->wrapTree.EvID.Value->Timestamp;
                first_time_valid = true;
            }

            timeInRun += t->wrapTree.Clock() / 1.0e6;

            auto evTime = (  timeInRun
                             + t->wrapTree.EvID.Value->Timestamp
                             - first_time );

            f(t, evTime);
        }
    }
}



treeLoader_t::treeLoader_t(const string& filename):
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

treeLoader_t::means_t treeLoader_t::getMeans() const
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

void taggEffTriple_t::initBkgFits()
{
    AvgBkgRates = timedData::getRatesVsTime({addressof(Bkg1),addressof(Bkg2)});
    double dummy(0);
    double tmax(0);
    AvgBkgRates->GetPoint(AvgBkgRates->GetN()-1,tmax,dummy);

    AvgBkgFit   = new TF1("avg","[0] + [1] * exp( - [2] * x)",0,tmax);
    AvgBkgFit->SetParameters(startparams.data());
    for (auto i: {0,1,2})
        AvgBkgFit->SetParLimits(i,0,upperlimit);

    AvgBkgRates->Fit(AvgBkgFit,"Q");
    for ( auto ch = 0u ; ch < Run.nchannels ; ++ch)
    {
        auto graph = timedData::getRatesVsTime({addressof(Bkg1),addressof(Bkg2)},ch);
        bkgFits.emplace_back(bkgFit_t(graph,ch));
        bkgFits.back().doFit(IntervalD(0,tmax),AvgBkgFit->GetParameter(2));
    }
}

taggEffTriple_t::taggEffTriple_t(const string& bkg1f, const string& runf, const string& bkg2f):
    Bkg1(bkg1f),
    Run(runf),
    Bkg2(bkg2f)
{
    LOG(INFO) << "Loading TaggEff measurement for:  "
              << "1st bkg: " << bkg1f << "  "
              << "run: "     << runf  << "  "
              << "2nd bkg: " << bkg2f;

    if ( !(Bkg1.setupName == Bkg2.setupName &&
           Bkg2.setupName == Run.setupName ) )
        throw runtime_error("Files in TaggEff-triple not from same Setup!");
    WrapTFileInput file(bkg2f);
    ant::TAntHeader* header;
    file.GetObject("AntHeader",header);
    if (!header)
        throw runtime_error("No Ant header in found!");
    startID = header->LastID;

    initBkgFits();

    avgRatesSub = new TGraph();
    avgRatesSub->SetMarkerStyle(kPlus);
    avgRatesSub->SetMarkerColor(kRed);

    avgRates = new TGraph();
    avgRates->SetMarkerStyle(kPlus);
    avgRates->SetMarkerColor(kGray);
}

unsigned taggEffTriple_t::sanityChecks(const treeLoader_t::means_t& bkg1, const treeLoader_t::means_t& run, const treeLoader_t::means_t& bkg2) const
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

const taggEff_t taggEffTriple_t::GetTaggEffSubtracted() const
{
    taggEff_t result(SetupName(),startID,Run.nchannels);

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

    vector<std_ext::RMS> TDCs(Run.nchannels);
    vector<std_ext::RMS> scalers(Run.nchannels);
    std_ext::RMS         L;

    auto time = 1.0 * Run.startTime - Bkg1.startTime;
    auto nentries = Run.wrapTree.Tree->GetEntries();

    for (auto entry = 0 ; entry < nentries ; ++entry)
    {
        Run.wrapTree.Tree->GetEntry(entry);
        time += Run.wrapTree.Clock() / 1.0e6;

        std_ext::RMS graphScaler;
        std_ext::RMS graphScalerSub;

        for (auto  ch = 0u ; ch < Run.nchannels ; ++ch)
        {
            TDCs.at(ch).Add(Run.wrapTree.TDCRates().at(ch));

            auto rate_sub = Run.wrapTree.TaggRates().at(ch) - bkgFits[ch](time);
            scalers.at(ch).Add( rate_sub );

            graphScaler.Add(Run.wrapTree.TaggRates().at(ch));
            graphScalerSub.Add(rate_sub);
        }
        FillGraph(avgRatesSub,time,graphScalerSub.GetMean());
        FillGraph(avgRates,time,graphScaler.GetMean());
        L.Add(Run.wrapTree.ExpLivetime);
    }

    for ( auto ch = 0u ; ch < Run.nchannels ; ++ch)
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
        result.BkgFitChi2.at(ch) = bkgFits.at(ch).Fit->GetChisquare();
    }


    return result;
}

const taggEff_t taggEffTriple_t::GetTaggEff() const
{
    taggEff_t result(SetupName(),startID,Run.nchannels);

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


    const auto m_bkg1 = Bkg1.getMeans();
    const auto m_run = Run.getMeans();
    const auto m_bkg2 = Bkg2.getMeans();
    if ( sanityChecks(m_bkg1,m_run,m_bkg2) > 0 )
        throw runtime_error("Tagg-Eff-triple didn't pass sanity checks, doublecheck input files!");

    auto nchannels = Bkg1.nchannels;
    result.TaggEffs.resize(nchannels,0);
    result.TaggEffErrors.resize(nchannels,0);
    result.BkgFitChi2.resize(nchannels,0);

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

taggEffTriple_t::bkgFit_t::bkgFit_t(TGraph* graph, const size_t channel): Graph(graph)
{
    string name = std_ext::formatter() << "channel" << channel;
    Graph->SetMarkerStyle(kPlus);
    Graph->SetMarkerColor(kBlue);
    Graph->SetTitle(name.c_str());
}

void taggEffTriple_t::bkgFit_t::doFit(const IntervalD& fitrange, const double lambda)
{
    string name = std_ext::formatter() << Graph->GetTitle() << "fit" ;
    string fitFkt = std_ext::formatter() << "[0] + [1] * exp( - " << lambda << " * x)";
    Fit   = new TF1(name.c_str(),fitFkt.c_str(),
                    fitrange.Start(),fitrange.Stop());
    for (auto i: {0,1})
    {
        Fit->SetParameter(i,startparams[i]);
        Fit->SetParLimits(i,0,upperlimit);
    }

    Graph->Fit(Fit,"Q");
}

double taggEffTriple_t::bkgFit_t::operator ()(const double time) const
{
    if ( Fit )
        return Fit->Eval(time);
    return std_ext::NaN;
}



TGraph* timedData::getRatesVsTime(const std::list<treeLoader_t*>& tContainers)
{
    auto graph2d = new TGraph();

    graph2d->SetMarkerStyle(kPlus);
    graph2d->GetXaxis()->SetTitle("time [s]");
    graph2d->GetYaxis()->SetTitle("avg. rate [Hz]");

    auto f = [graph2d] (treeLoader_t* t, double evTime) {
        std_ext::RMS rmsRate;
        for ( const auto& tr: t->wrapTree.TaggRates())
            rmsRate.Add(tr);
        FillGraph(graph2d,evTime, rmsRate.GetMean());
    };

    runOverContainer(tContainers, f);
    return graph2d;
}

TGraph* timedData::getRatesVsTime(const std::list<treeLoader_t*>& tContainers, const size_t channel)
{
    auto graph2d = new TGraph();

    graph2d->SetMarkerStyle(kPlus);
    graph2d->GetXaxis()->SetTitle("time [s]");
    graph2d->GetYaxis()->SetTitle("avg. rate [Hz]");

    auto f = [channel, graph2d] (treeLoader_t* t, double evTime) {
        FillGraph(graph2d,evTime,t->wrapTree.TaggRates().at(channel));
    };

    runOverContainer(tContainers, f);
    return graph2d;
}

TGraph* timedData::getLtVsTime(const std::list<treeLoader_t*>& tContainers)
{
    auto graph = new TGraph();

    graph->SetMarkerStyle(kPlus);
    graph->GetXaxis()->SetTitle("time [s]");
    graph->GetYaxis()->SetTitle("livetime");

    auto f = [graph] (treeLoader_t* t, double evTime)
    {
        FillGraph(graph,evTime,t->wrapTree.ExpLivetime);
    };
    runOverContainer(tContainers,f);
    return graph;
}
