#include "base/CmdLine.h"
#include "base/Logger.h"

#include "base/WrapTFile.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "TTree.h"

#include <limits>

using namespace ant;
using namespace std;


int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("compare_ant_goat", ' ', "0.1");

    auto cmd_antinput = cmd.add<TCLAP::ValueArg<string>>("","ant","treeEvents from Ant",true,"","rootfile");
    auto cmd_goatinput = cmd.add<TCLAP::ValueArg<string>>("","goat","treeEvents from Goat",true,"","rootfile");

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
    long long nEvents = 0;

    while(entryAnt < treeAnt->GetEntries() && entryGoat < treeGoat->GetEntries()) {

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

        if(antRecon.TaggerHits.size() != goatRecon.TaggerHits.size()) {
            LOG(ERROR) << "TaggerHits size mismatch at ant=" << entryAnt << " goat=" << entryGoat;
            break;
        }

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
                if(!isfinite(hit.Energy) || !isfinite(hit.Time))
                    continue;
                cout << "   Ch=" << it_hit.first
                     <<  " Energy=" << hit.Energy << " Time=" << hit.Time << endl;
            }
        };

        cout << "<<<< Ant Event=" << nEvents << endl;
        print_dethits(antRecon, Detector_t::Type_t::CB);
        print_dethits(antRecon, Detector_t::Type_t::PID);

        cout << endl;

        cout << "<<<< Goat Event=" << nEvents << endl;
        print_dethits(goatRecon, Detector_t::Type_t::CB);
        print_dethits(goatRecon, Detector_t::Type_t::PID);

        cout << endl;

        auto compare_clusters = [] (
                                const TEventData& antRecon,
                                const TEventData& goatRecon,
                                Detector_t::Type_t type,
                                double threshold
                                ) {

            // compare clusters Ant/Goat
            struct cluster_t {
                TClusterList::const_iterator Ant;
                TClusterList::const_iterator Goat;
            };

            map<unsigned, cluster_t> clusters;
            double antEnergy = 0;
            for(auto it_cl = antRecon.Clusters.begin(); it_cl; ++it_cl)
                if(it_cl->DetectorType == type && it_cl->Energy>threshold) {
                    antEnergy += it_cl->Energy;
                    clusters[it_cl->CentralElement].Ant = it_cl;
                }

            double goatEnergy = 0;
            for(auto it_cl = goatRecon.Clusters.begin(); it_cl; ++it_cl)
                if(it_cl->DetectorType == type && it_cl->Energy>threshold) {
                    goatEnergy += it_cl->Energy;
                    clusters[it_cl->CentralElement].Goat = it_cl;
                }

            cout << "  EnergySum: " << antEnergy << "/" << goatEnergy << endl;
            for(const auto& it_cl : clusters) {
                auto stringify_cl = [] (const TClusterList::const_iterator& cl) -> string {
                    if(!cl)
                        return "";
                    return std_ext::formatter() << "(E=" << cl->Energy
                                                << (cl->HasFlag(TCluster::Flags_t::Split) ? ",S" : "")
                                                << ")";
                };
                const cluster_t& cl = it_cl.second;

                cout << "  Ch=" << it_cl.first
                     << " Ant=" << stringify_cl(cl.Ant)
                     << " Goat=" << stringify_cl(cl.Goat);
                if(cl.Ant && cl.Goat)
                    cout << " Angle=" << cl.Ant->Position.Angle(cl.Goat->Position);
                cout << endl;
            }
        };

        cout << "<<<< CB Clusters Event=" << nEvents << endl;
        compare_clusters(antRecon, goatRecon, Detector_t::Type_t::CB, 25);

        cout << endl;

        cout << "<<<< PID Clusters Event=" << nEvents << endl;
        compare_clusters(antRecon, goatRecon, Detector_t::Type_t::PID, 0);

        cout << endl;

        // compare candidates Ant/Goat
        struct candidate_t {
            TCandidatePtr Ant;
            TCandidatePtr Goat;
        };

        map<unsigned, candidate_t> candidates;

        for(const TCandidatePtr& c : antRecon.Candidates)
            if(c->Detector & Detector_t::Any_t::CB_Apparatus)
                candidates[c->FindCaloCluster()->CentralElement].Ant = c;

        for(const TCandidatePtr& c : goatRecon.Candidates)
            if(c->Detector & Detector_t::Any_t::CB_Apparatus)
                candidates[c->FindCaloCluster()->CentralElement].Goat = c;

        cout << "<<<< Candidates Event=" << nEvents << endl;
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

        cout << endl << endl;

        nEvents++;
        if(nEvents==50)
            break;
        entryAnt++;
        entryGoat++;
    }

    LOG(INFO) << "Compared " << nEvents << " events";
}