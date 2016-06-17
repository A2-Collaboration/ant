#include "base/CmdLine.h"
#include "base/Logger.h"

#include "base/WrapTFile.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "TTree.h"
#include "TRint.h"
#include "TH1D.h"
#include "TH2D.h"


#include "analysis/plot/HistogramFactories.h"
#include "analysis/plot/root_draw.h"

#include <limits>

using namespace ant;
using namespace std;
using namespace ant::analysis;
volatile bool interrupt = false;


int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        cout << ">>> Interrupted" << endl;
        interrupt = true;
    });

    TCLAP::CmdLine cmd("compare_ant_goat", ' ', "0.1");

    auto cmd_antinput = cmd.add<TCLAP::ValueArg<string>>("","ant","treeEvents from Ant",true,"","rootfile");
    auto cmd_goatinput = cmd.add<TCLAP::ValueArg<string>>("","goat","treeEvents from Goat",true,"","rootfile");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_dump = cmd.add<TCLAP::MultiSwitchArg>("","dump","Dump events to console",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");


    cmd.parse(argc, argv);

    WrapTFileInput antinput(cmd_antinput->getValue());
    WrapTFileInput goatinput(cmd_goatinput->getValue());

    TTree* treeAnt;
    if(!antinput.GetObject("treeEvents", treeAnt)) {
        LOG(ERROR) << "Cannot find treeEvents in " << cmd_antinput->getValue();
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "Ant treeEvent Entries " << treeAnt->GetEntries();

    TTree* treeGoat;
    if(!goatinput.GetObject("treeEvents", treeGoat)) {
        LOG(ERROR) << "Cannot find treeEvents in " << cmd_goatinput->getValue();
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "Goat treeEvent Entries " << treeGoat->GetEntries();


    TEvent* eventAnt = new TEvent();
    TEvent* eventGoat = new TEvent();

    treeAnt->SetBranchAddress("data",addressof(eventAnt));
    treeGoat->SetBranchAddress("data",addressof(eventGoat));

    long long entryAnt = 0;
    long long entryGoat = 0;

    bool synced = false;
    const bool dump = cmd_dump->isSet();
    long long maxevents = cmd_maxevents->isSet()
            ? cmd_maxevents->getValue().back()
            :  numeric_limits<long long>::max();
    long long nEvents = 0;

    HistogramFactory HistFac("compare_ant_goat");

    BinSettings bins_CaloEnergy(200,0,600);
    auto h_CaloEnergy_CB = HistFac.makeTH2D("h_CaloEnergy CB","Goat","Ant",bins_CaloEnergy,bins_CaloEnergy,"h_CaloEnergy_CB");
    auto h_CaloEnergy_TAPS = HistFac.makeTH2D("h_CaloEnergy TAPS","Goat","Ant",bins_CaloEnergy,bins_CaloEnergy,"h_CaloEnergy_TAPS");

    BinSettings bins_VetoEnergy(30,0,10);
    auto h_VetoEnergy_CB = HistFac.makeTH2D("h_VetoEnergy CB","Goat","Ant",bins_VetoEnergy,bins_VetoEnergy,"h_VetoEnergy_CB");
    auto h_VetoEnergy_TAPS = HistFac.makeTH2D("h_VetoEnergy TAPS","Goat","Ant",bins_VetoEnergy,bins_VetoEnergy,"h_VetoEnergy_TAPS");


    while(entryAnt < treeAnt->GetEntries() && entryGoat < treeGoat->GetEntries()) {

        if(interrupt)
            break;

        treeAnt->GetEntry(entryAnt);
        treeGoat->GetEntry(entryGoat);

        if(!synced) {
            const auto antEventID = eventAnt->Reconstructed().Trigger.DAQEventID;
            const auto goatEventID = eventGoat->Reconstructed().Trigger.DAQEventID;
            if(antEventID < goatEventID) {
                entryAnt++;
                continue;
            }
            else if(antEventID > goatEventID) {
                entryGoat++;
                continue;
            }
            synced = true;
            LOG(INFO) << "Synced at ant=" << entryAnt << " goat=" << entryGoat;
        }

        const TEventData& antRecon = eventAnt->Reconstructed();
        const TEventData& goatRecon = eventGoat->Reconstructed();

//        if(antRecon.TaggerHits.size() != goatRecon.TaggerHits.size()) {
//            LOG(ERROR) << "TaggerHits size mismatch at ant=" << entryAnt << " goat=" << entryGoat;
//            break;
//        }

        // only events with 2 clusters in CB
//        auto count_cbclusters = [] (const TEventData& recon) {
//            unsigned nCBClusters = 0;
//            for(const auto& cl : recon.Clusters) {
//                if(cl->DetectorType == Detector_t::Type_t::CB)
//                    nCBClusters++;
//            }
//            return nCBClusters;
//        };

//        if(count_cbclusters(antRecon) != 2 && count_cbclusters(goatRecon) != 2) {
//            entryAnt++;
//            entryGoat++;
//            continue;
//        }

        if(dump) {

            const auto print_dethits = [] (const TEventData& recon, Detector_t::Type_t type) {
                // match hits by channel again
                struct hit_t {
                    double Energy = numeric_limits<double>::quiet_NaN();
                    double Time = numeric_limits<double>::quiet_NaN();
                };
                map<unsigned, hit_t> hits;
                for(const TDetectorReadHit& dethit : recon.DetectorReadHits) {
                    if(dethit.Values.empty())
                        continue;
                    if(dethit.DetectorType  != type)
                        continue;
                    hit_t& hit = hits[dethit.Channel];
                    if(dethit.ChannelType == Channel_t::Type_t::Integral)
                        hit.Energy = dethit.Values.back();
                    if(dethit.ChannelType == Channel_t::Type_t::Timing)
                        hit.Time = dethit.Values.back();
                }

                cout << ">> DetectorHits n=" << hits.size() << " " << Detector_t::ToString(type) << endl;
                for(const auto& it_hit : hits) {
                    const hit_t& hit = it_hit.second;
                    //                if(!isfinite(hit.Energy) || !isfinite(hit.Time))
                    //                    continue;
                    cout << "   Ch=" << it_hit.first
                         <<  " Energy=" << hit.Energy << " Time=" << hit.Time << endl;
                }
            };


            cout << "<<<< Ant Event=" << nEvents << endl;
            //        print_dethits(antRecon, Detector_t::Type_t::CB);
            print_dethits(antRecon, Detector_t::Type_t::PID);

            cout << endl;

            cout << "<<<< Goat Event=" << nEvents << endl;
            //        print_dethits(goatRecon, Detector_t::Type_t::CB);
            print_dethits(goatRecon, Detector_t::Type_t::PID);

            cout << endl;

        }

//        auto compare_clusters = [] (
//                                const TEventData& antRecon,
//                                const TEventData& goatRecon,
//                                Detector_t::Type_t type,
//                                double threshold
//                                ) {

//            // compare clusters Ant/Goat
//            struct cluster_t {
//                TClusterPtr Ant;
//                TClusterPtr Goat;
//            };

//            map<unsigned, cluster_t> clusters;
//            double antEnergy = 0;
//            for(auto it_cl = antRecon.Clusters.begin(); it_cl != antRecon.Clusters.end(); ++it_cl)
//                if(it_cl->DetectorType == type && it_cl->Energy>threshold) {
//                    antEnergy += it_cl->Energy;
//                    clusters[it_cl->CentralElement].Ant = it_cl.get_ptr();
//                }

//            double goatEnergy = 0;
//            for(auto it_cl = goatRecon.Clusters.begin(); it_cl != goatRecon.Clusters.end(); ++it_cl)
//                if(it_cl->DetectorType == type && it_cl->Energy>threshold) {
//                    goatEnergy += it_cl->Energy;
//                    clusters[it_cl->CentralElement].Goat = it_cl.get_ptr();
//                }

//            cout << "  EnergySum: " << antEnergy << "/" << goatEnergy << endl;
//            for(const auto& it_cl : clusters) {
//                auto stringify_cl = [] (const TClusterPtr& cl) -> string {
//                    if(!cl)
//                        return "";
//                    return std_ext::formatter() << "(E=" << cl->Energy
//                                                << (cl->HasFlag(TCluster::Flags_t::Split) ? ",S" : "")
//                                                << ")";
//                };
//                const cluster_t& cl = it_cl.second;

//                cout << "  Ch=" << it_cl.first
//                     << " Ant=" << stringify_cl(cl.Ant)
//                     << " Goat=" << stringify_cl(cl.Goat);
//                if(cl.Ant && cl.Goat)
//                    cout << " Angle=" << cl.Ant->Position.Angle(cl.Goat->Position);
//                cout << endl;
//            }
//        };

//        cout << "<<<< CB Clusters Event=" << nEvents << endl;
//        compare_clusters(antRecon, goatRecon, Detector_t::Type_t::CB, 25);

//        cout << endl;

//        cout << "<<<< PID Clusters Event=" << nEvents << endl;
//        compare_clusters(antRecon, goatRecon, Detector_t::Type_t::PID, 0);

//        cout << endl;

        // compare candidates Ant/Goat
        struct candidate_t {
            TCandidatePtr Ant;
            TCandidatePtr Goat;
        };

        using mapped_candidates_t = map<unsigned, candidate_t>;
        mapped_candidates_t candidates_cb;
        mapped_candidates_t candidates_taps;


        for(const auto& c : antRecon.Candidates.get_iter()) {
            if(c->Detector & Detector_t::Any_t::CB_Apparatus)
                candidates_cb[c->FindCaloCluster()->CentralElement].Ant = c;
            else if(c->Detector & Detector_t::Any_t::TAPS_Apparatus)
                candidates_taps[c->FindCaloCluster()->CentralElement].Ant = c;
        }

        for(const auto& c : goatRecon.Candidates.get_iter()) {
            if(c->Detector & Detector_t::Any_t::CB_Apparatus)
                candidates_cb[c->FindCaloCluster()->CentralElement].Goat = c;
            else if(c->Detector & Detector_t::Any_t::TAPS_Apparatus)
                candidates_taps[c->FindCaloCluster()->CentralElement].Goat = c;
        }


        if(dump) {
            auto dump_candidates = [] (const mapped_candidates_t& candidates) {
                for(const auto& it_c : candidates) {
                    auto stringify_c = [] (const TCandidatePtr& c) -> string {
                        if(!c)
                            return "";
                        return std_ext::formatter()
                                << "(E=" << c->CaloEnergy
                                << ",Cl=" << c->ClusterSize
                                << ",vE=" << c->VetoEnergy
                                << ",t=" << c->Time
                                << ")";
                    };
                    const candidate_t& c = it_c.second;
                    cout << "  Ch=" << it_c.first
                         << " Ant=" << stringify_c(c.Ant)
                         << " Goat=" << stringify_c(c.Goat)
                         << endl;
                }
            };


            cout << "<<<< Candidates CB Event=" << nEvents << endl;
            dump_candidates(candidates_cb);
            cout << endl << endl;

            cout << "<<<< Candidates TAPS Event=" << nEvents << endl;
            dump_candidates(candidates_taps);
            cout << endl << endl;
        }


        for(const auto& it_c : candidates_cb) {
            const candidate_t& c = it_c.second;
            if(!c.Ant || !c.Goat)
                continue;
            h_CaloEnergy_CB->Fill(c.Ant->CaloEnergy, c.Goat->CaloEnergy);
            h_VetoEnergy_CB->Fill(c.Ant->VetoEnergy, c.Goat->VetoEnergy);

        }

        for(const auto& it_c : candidates_taps) {
            const candidate_t& c = it_c.second;
            if(!c.Ant || !c.Goat)
                continue;
            h_CaloEnergy_TAPS->Fill(c.Ant->CaloEnergy, c.Goat->CaloEnergy);
            h_VetoEnergy_TAPS->Fill(c.Ant->VetoEnergy, c.Goat->VetoEnergy);
        }

        nEvents++;
        if(nEvents==maxevents)
            break;
        entryAnt++;
        entryGoat++;
    }

    LOG(INFO) << "Compared " << nEvents << " events";

    if(!cmd_batchmode->isSet()) {

        argc=0; // prevent TRint to parse any cmdline
        TRint app("Ant",&argc,argv,nullptr,0,true);

        canvas("compare_ant_goat") << drawoption("colz")
                                   << h_CaloEnergy_CB
                                   << h_CaloEnergy_TAPS
                                   << h_VetoEnergy_CB
                                   << h_VetoEnergy_TAPS
                                   << endc;

        app.Run(kTRUE); // really important to return...

    }
}