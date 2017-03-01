#include "TAPSEnergy_Check.h"

#include "analysis/plot/HistogramFactory.h"
#include "expconfig/ExpConfig.h"
#include "analysis/plot/root_draw.h"
#include "base/interval.h"
#include "base/std_ext/math.h"
#include "tree/TParticle.h"

#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"


#include <iostream>
#include <sstream>

using namespace std;
using namespace ant;
using namespace ant::analysis;



void TAPSEnergy_Check::Analyse(TFile* file)
{
    auto h = dynamic_cast<TH3D*>(file->Get("TAPS_Energy/ggIM_mult"));

    canvas c;
    for(int m=0;m<=10;m++) {
        h->GetZaxis()->SetRangeUser(m,m+1);
        stringstream ss_name;
        ss_name << "Mult_" << m << "_yx";
        c << drawoption("colz") << h->Project3D(ss_name.str().c_str());

    }
    c << dynamic_cast<TH2D*>(file->Get("TAPS_Energy/ggIM"));
    c << endc;
}

void TAPSEnergy_Check::AnalyseTree(TFile* file)
{
    TTree* tree = dynamic_cast<TTree*>(file->Get("TAPS_Energy/cands_tree"));
    struct tree_data_t {
        vector<double>* Ek = nullptr;
        vector<double>* Theta = nullptr;
        vector<double>* Phi = nullptr;
        vector<double>* VetoE = nullptr;
        vector<double>* Time = nullptr;
        vector<unsigned>* Channel = nullptr;
        tree_data_t(const std::string& prefix, TTree* tree) {
            tree->SetBranchAddress((prefix+"_Ek").c_str(), &Ek);
            tree->SetBranchAddress((prefix+"_Theta").c_str(), &Theta);
            tree->SetBranchAddress((prefix+"_Phi").c_str(), &Phi);
            tree->SetBranchAddress((prefix+"_VetoE").c_str(), &VetoE);
            tree->SetBranchAddress((prefix+"_Time").c_str(), &Time);
            tree->SetBranchAddress((prefix+"_Channel").c_str(), &Channel);

        }
    };
    tree_data_t CB("CB",tree);
    tree_data_t TAPS("TAPS",tree);

    struct hist_t {
        TH1D* IM;
        TH1D* nTAPS;
        TH1D* nCB;
        hist_t(unsigned ch) {
            string titlestr = std_ext::formatter() << "TAPSEnergy_Ch" << ch;
            HistogramFactory HistFac(titlestr);
            HistFac.SetTitlePrefix(titlestr);
            IM = HistFac.makeTH1D("IM","IM / MeV","#",BinSettings(400,0,500),"IM");
            nTAPS = HistFac.makeTH1D("nTAPS","Photons in TAPS","#",BinSettings(10),"nTAPS");
            nCB = HistFac.makeTH1D("nCB","Photons in CB","#",BinSettings(10),"nCB");
        }
    };

    map<unsigned, hist_t> hists;

    for(Long64_t entry=0;entry<tree->GetEntries();entry++)  {
        tree->GetEntry(entry);

        std::list<std::pair<unsigned, TParticle>> taps_photons;
        std::list<TParticle> cb_photons;

        for(unsigned i=0;i<TAPS.Channel->size();i++) {
            // cut on TAPS properties
            if(fabs((*TAPS.Time)[i]) > 5)
                continue;
            if((*TAPS.VetoE)[i] > 0.5)
                continue;

            // build the TAPS photon but remember the channel
            TParticle photon(ParticleTypeDatabase::Photon, (*TAPS.Ek)[i], (*TAPS.Theta)[i], (*TAPS.Phi)[i]);
            taps_photons.emplace_back((*TAPS.Channel)[i], move(photon));
        }

        for(unsigned j=0;j<CB.Channel->size();j++) {
            if((*CB.VetoE)[j] > 0.5)
                continue;

            cb_photons.emplace_back(ParticleTypeDatabase::Photon, (*CB.Ek)[j], (*CB.Theta)[j], (*CB.Phi)[j]);
        }

        if(cb_photons.size() < 2 || cb_photons.size() > 6)
            continue;

        for(const auto& p_taps_photon : taps_photons) {
            // find the histogram by TAPS channel
            unsigned ch = p_taps_photon.first;
            auto it_h = hists.find(ch);
            if(it_h == hists.end()) {
               auto p = hists.emplace(make_pair(ch, hist_t(ch)));
               if(p.second)
                   it_h = p.first;
            }
            hist_t& h = it_h->second;

            h.nTAPS->Fill(taps_photons.size());
            h.nCB->Fill(cb_photons.size());

            const auto& taps_photon = p_taps_photon.second;
            for(const auto& cb_photon : cb_photons) {
                const auto pi0 = taps_photon + cb_photon;
                h.IM->Fill(pi0.M());
            }
        }



    }


    canvas IM("IM");
    canvas nTAPS("nTAPS");
    nTAPS << padoption::enable(padoption::LogY);
    canvas nCB("nCB");
    nCB << padoption::enable(padoption::LogY);
    for(auto it_h : hists) {
        IM << it_h.second.IM;
        nTAPS << it_h.second.nTAPS;
        nCB << it_h.second.nCB;
    }

    nTAPS << endc;
    nCB << endc;
    IM << endc;
}

