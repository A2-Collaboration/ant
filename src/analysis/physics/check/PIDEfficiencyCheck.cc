#include "PIDEfficiencyCheck.h"

#include "physics/Plotter.h"
#include "utils/uncertainties/Interpolated.h"
#include "plot/CutTree.h"

#include "expconfig/ExpConfig.h"

#include "base/ParticleTypeTree.h"
#include "base/std_ext/misc.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::plot;

PIDEfficiencyCheck::PIDEfficiencyCheck(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    fitter(utils::UncertaintyModels::Interpolated::makeAndLoad(),
           true // enable fit z vertex
           )
{
    promptrandom.AddPromptRange({-2.5,2.5});
    promptrandom.AddRandomRange({-35,-2.5});
    promptrandom.AddRandomRange({ 2.5,35});
    fitter.SetZVertexSigma(3.0);
    t.CreateBranches(HistFac.makeTTree("t"));

    // set some iteration filter to speed up fitting
//    auto pi0 = fitter.GetTreeNode(ParticleTypeDatabase::Pi0);
//    fitter.SetIterationFilter(
//                [pi0] () {
//        return ParticleTypeDatabase::Pi0.GetWindow(70).Contains(pi0->Get().LVSum.M());
//    });

    h_Steps = HistFac.makeTH1D("Steps","","",BinSettings(10),"h_Steps");
}

void PIDEfficiencyCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    h_Steps->Fill("Seen",1.0);

    const auto& data = event.Reconstructed();

    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h_Steps->Fill("TagHits total", 1.0);
    }

    if(data.Candidates.size() != 3)
        return;
    h_Steps->Fill("nCands==3",1.0);

    // fill some PID info
    t.PID_Ch().clear();
    t.PID_Phi().clear();
    t.PID_E().clear();
    t.PID_Time().clear();
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType != Detector_t::Type_t::PID)
            continue;
        t.PID_Ch().push_back(cl.CentralElement);
        t.PID_Phi().push_back(std_ext::radian_to_degree(cl.Position.Phi()));
        t.PID_E().push_back(cl.Energy);
        t.PID_Time().push_back(cl.Time);
    }


    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h_Steps->Fill("TagHits",1.0);
        h_Steps->Fill("TagHits prompt",0.0);
        if(promptrandom.State() == PromptRandom::Case::Prompt)
            h_Steps->Fill("TagHits prompt",1.0);


        t.TaggCh = taggerhit.Channel;
        t.TaggE  = taggerhit.PhotonEnergy;
        t.TaggW  = promptrandom.FillWeight();
        t.TaggT  = taggerhit.Time;
        t.CBAvgTime = triggersimu.GetRefTiming();

        t.FitProb = std_ext::NaN;

        for(auto cand_proton : data.Candidates.get_iter()) {
            auto proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);
            TParticleList photons;
            LorentzVec photon_sum;
            for(auto cand_photon : data.Candidates.get_iter()) {
                if(cand_photon == cand_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_photon));
                photon_sum += *photons.back();
            }
            if(photon_sum.M()>250)
                continue;

            APLCON::Result_t fitresult = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

            if(fitresult.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(t.FitProb, fitresult.Probability))
                continue;
            auto fitted_proton = fitter.GetFittedProton();
            auto fitted_photons = fitter.GetFittedPhotons();

            t.Proton_Ek = fitted_proton->Ek();
            t.Proton_Phi = std_ext::radian_to_degree(fitted_proton->Phi());
            t.Proton_Theta = std_ext::radian_to_degree(fitted_proton->Theta());
            t.Proton_VetoE = fitted_proton->Candidate->VetoEnergy;

            t.Proton_inCB = (bool)(proton->Candidate->Detector & Detector_t::Type_t::CB);

            t.Proton_MinPIDPhi = std_ext::NaN;
            for(unsigned i=0;i<t.PID_Phi().size();i++) {
                const auto diff = t.Proton_Phi - t.PID_Phi().at(i);
                if(isfinite(t.Proton_MinPIDPhi) && abs(t.Proton_MinPIDPhi) < abs(diff))
                    continue;
                t.Proton_MinPIDPhi = diff;
                t.Proton_MinPIDCh = t.PID_Ch().at(i);
            }

            t.nPhotonsCB = (bool)(photons.back()->Candidate->Detector & Detector_t::Type_t::CB) +
                           (bool)(photons.front()->Candidate->Detector & Detector_t::Type_t::CB);

            t.IM_2g = photon_sum.M();
            t.IM_2g_fitted = (*fitted_photons.front() + *fitted_photons.back()).M();

        }

        if(isfinite(t.FitProb)) {
            h_Steps->Fill("Fit Ok",1.0);
            h_Steps->Fill("Proton CB", t.Proton_inCB);
            h_Steps->Fill(">0 Photon CB", t.nPhotonsCB>0);
        }


        // require reasonable fit and proton in CB
        if(t.FitProb>0.01) {
            h_Steps->Fill("Fills",1.0);
            if(t.Proton_VetoE>0)
                h_Steps->Fill("ProtonVetoE>0",1.0);
            t.Tree->Fill();
        }
    }
}

void PIDEfficiencyCheck::ShowResult()
{
    canvas(GetName())
            << h_Steps
            << drawoption("colz")
            << TTree_drawable(t.Tree, "Proton_MinPIDCh:Proton_MinPIDPhi >> (100,-70,70,24,0,24)","")
            << TTree_drawable(t.Tree, "Proton_VetoE:Proton_Ek","Proton_VetoE>0")
            << TTree_drawable(t.Tree, "TaggT-CBAvgTime","")
            << TTree_drawable(t.Tree, "Proton_inCB","")
            << endc;
}

struct Hist_t {
    using Tree_t = PIDEfficiencyCheck::Tree_t;
    struct Fill_t {
        Fill_t(const Tree_t& t) : Tree(t) {}
        const Tree_t& Tree;
        double Weight() const {
            return Tree.TaggW;
        }
    };

    TH1D* h_KinFitProb;
    TH1D* h_TaggT;
    TH1D* h_IM_2g;
    TH1D* h_IM_2g_fitted;
    TH1D* h_IM_2g_fitted_sub;
    TH1D* h_nPhotonsCB;
    TH2D* h_ProtonThetaEk;
    TH2D* h_ProtonMinPIDPhi;
    TH1D* h_ProtonMinPIDCh;
    TH1D* h_ProtonPhi;
    TH2D* h_Banana;


    Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t) {
        BinSettings bins_IM(100,50,250);

        auto pid = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);

        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",BinSettings(100,0,1),"h_KinFitProb");
        h_TaggT = HistFac.makeTH1D("TaggT-CBAvgTime","t / ns","",BinSettings(200,-40,40),"h_TaggT");
        h_IM_2g = HistFac.makeTH1D("IM_2g","IM / MeV","",bins_IM,"h_IM_2g");
        h_IM_2g_fitted = HistFac.makeTH1D("IM_2g_fitted","IM / MeV","",bins_IM,"h_IM_2g_fitted");
        h_IM_2g_fitted_sub = HistFac.makeTH1D("IM_2g_fitted_sub","IM / MeV","",bins_IM,"h_IM_2g_fitted_sub");
        h_nPhotonsCB = HistFac.makeTH1D("nPhotonsCB","","",BinSettings(3),"h_nPhotonsCB");
        h_ProtonThetaEk = HistFac.makeTH2D("Proton #theta E_k","Ek / MeV","#theta / #circ",
                                           BinSettings(100,0,1000),BinSettings(40,20,80),
                                           "h_ProtonThetaEk");
        h_ProtonMinPIDPhi = HistFac.makeTH2D("Proton MinPID","#Delta#phi / #circ","PID Channel",
                                             BinSettings(200,-200,200),BinSettings(pid->GetNChannels()),
                                             "h_ProtonMinPIDPhi");
        h_ProtonMinPIDCh = HistFac.makeTH1D("Proton MinPIDCh","PID Channel","",
                                            BinSettings(pid->GetNChannels()),
                                            "h_ProtonMinPIDCh");
        h_ProtonPhi = HistFac.makeTH1D("Proton Phi","#phi / #circ","",BinSettings(pid->GetNChannels(), -180, 180),"h_ProtonPhi");
        h_Banana = HistFac.makeTH2D("Banana","E_k / MeV","VetoE / MeV",
                                    BinSettings(100,0,1000),BinSettings(100,0,10),
                                    "h_Banana");
    }

    void Fill(const Fill_t& f) const {
        h_KinFitProb->Fill(f.Tree.FitProb);
        h_TaggT->Fill(f.Tree.TaggT-f.Tree.CBAvgTime);
        h_IM_2g->Fill(f.Tree.IM_2g);
        h_IM_2g_fitted->Fill(f.Tree.IM_2g_fitted);
        h_IM_2g_fitted_sub->Fill(f.Tree.IM_2g_fitted, f.Weight());
        h_nPhotonsCB->Fill(f.Tree.nPhotonsCB);
        h_ProtonThetaEk->Fill(f.Tree.Proton_Ek, f.Tree.Proton_Theta);
        h_ProtonMinPIDPhi->Fill(f.Tree.Proton_MinPIDPhi,f.Tree.Proton_MinPIDCh);
        h_ProtonMinPIDCh->Fill(f.Tree.Proton_MinPIDCh);
        h_ProtonPhi->Fill(f.Tree.Proton_Phi);
        h_Banana->Fill(f.Tree.Proton_Ek, f.Tree.Proton_VetoE);
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        using i_t = interval<double>;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"prompt", [] (const Fill_t& f) { return i_t(-2.5, 2.5).Contains(f.Tree.TaggT-f.Tree.CBAvgTime); }},
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"100<IM<180", [] (const Fill_t& f) { return i_t(100, 180).Contains(f.Tree.IM_2g_fitted); }},
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"ProtonTheta>40", [] (const Fill_t& f) { return f.Tree.Proton_Theta>40; } },
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"-25<Phi<25", [] (const Fill_t& f) { return i_t(-25,25).Contains(f.Tree.Proton_MinPIDPhi); } },
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        return cuts;
    }
};

struct PIDEfficiencyCheck_plot : Plotter {

    Hist_t::Tree_t tree;

    cuttree::Tree_t<Hist_t> mycuttree;

    PIDEfficiencyCheck_plot(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        Plotter(name, input, opts)
    {
        if(!input.GetObject("PIDEfficiencyCheck/t", tree.Tree))
            throw Exception("Cannot find tree PIDEfficiencyCheck/t");
        tree.LinkBranches();

        mycuttree = cuttree::Make<Hist_t>(HistFac);
    }

    virtual long long GetNumEntries() const override
    {
        return tree.Tree->GetEntries();
    }

    virtual void ProcessEntry(const long long entry) override
    {
        tree.Tree->GetEntry(entry);
        cuttree::Fill<Hist_t>(mycuttree, tree);
    }

};

AUTO_REGISTER_PHYSICS(PIDEfficiencyCheck)
AUTO_REGISTER_PLOTTER(PIDEfficiencyCheck_plot)
