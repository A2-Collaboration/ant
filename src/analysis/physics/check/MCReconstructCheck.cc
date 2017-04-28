#include "MCReconstructCheck.h"

#include "utils/ParticleTools.h"
#include "root-addons/analysis_codes/hstack.h"
#include "root-addons/cbtaps_display/TH2TAPS.h"
#include "plot/HistStyle.h"
#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/TAPS.h"
#include "base/Detector_t.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TH1D.h"
#include "TTree.h"

#include <cmath>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCReconstructCheck::MCReconstructCheck(const std::string& name, OptionsPtr opts):
    Physics(name,opts),
    cb_group(HistFac, "CB", histgroup::detectortype::CB),
    taps_group(HistFac,"TAPS",histgroup::detectortype::TAPS),
    all_group(HistFac,"All",histgroup::detectortype::All),
    tapsveto(HistFac),
    mult1_only(opts->Get<bool>("Mult1Only",false)),
    show_cb_only(opts->Get<bool>("ShowCBonly",false)),
    fill_tree(opts->Get<bool>("FillTree",false))
{
    if(mult1_only)
        LOG(INFO) << "Using multiplicity == 1 events only";

    if(fill_tree)
        t.CreateBranches(HistFac.makeTTree("tree"));

    timesmear.smearing_enabled = true;
}

Detector_t::Any_t GetCommonDetector(const TCandidateList& cands) {
    Detector_t::Any_t common_detetctor = Detector_t::Any_t::None;
    for(const auto& cand : cands) {
        common_detetctor |= cand.Detector;
    }
    return common_detetctor;
}

void MCReconstructCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);

    if(mctrue_particles.GetAll().size() == 1) {

        if(mult1_only && event.Reconstructed().Candidates.size() > 1)
            return;

        const auto& mctrue_particle = mctrue_particles.GetAll().at(0);

        const auto common_detector = GetCommonDetector(event.Reconstructed().Candidates);

        if(common_detector & Detector_t::Any_t::CB_Apparatus) {
            cb_group.Fill(mctrue_particle, event.Reconstructed().Candidates, event.Reconstructed().Clusters);
        }

        if(common_detector & Detector_t::Any_t::TAPS_Apparatus) {
            taps_group.Fill(mctrue_particle, event.Reconstructed().Candidates, event.Reconstructed().Clusters);
        }

        all_group.Fill(mctrue_particle, event.Reconstructed().Candidates, event.Reconstructed().Clusters);

        tapsveto.Fill(event.Reconstructed().Candidates, event.Reconstructed().Clusters);

        t.mult = unsigned(event.Reconstructed().Candidates.size());

        t.tE     = mctrue_particle->Ek();
        t.tTheta = std_ext::radian_to_degree(mctrue_particle->Theta());
        t.tPhi   = std_ext::radian_to_degree(mctrue_particle->Phi());

        for(const auto& c : event.Reconstructed().Candidates) {
            t.rE     = c.CaloEnergy;
            t.rTheta = std_ext::radian_to_degree(c.Theta);
            t.rPhi   = std_ext::radian_to_degree(c.Phi);
            t.rVeto  = c.VetoEnergy;
            t.rTime  = timesmear.GetTime(c);
            t.rSize  = c.ClusterSize;
            if(c.Detector & Detector_t::Any_t::CB_Apparatus)
                t.Cal    = 1;
            else if(c.Detector & Detector_t::Any_t::TAPS_Apparatus)
                t.Cal    = 2;
            else
                t.Cal    = 0;

            if(t)
                t.Tree->Fill();
        }



    }
}


void MCReconstructCheck::Finish()
{
    cb_group.Finish();
    taps_group.Finish();
    all_group.Finish();
    tapsveto.Finish();
    LOG(INFO) << "MCReconstructCheck Finish";
}

void MCReconstructCheck::ShowResult()
{
    cb_group.ShowResult();
    if(show_cb_only)
        return;
    taps_group.ShowResult();
    all_group.ShowResult();
    tapsveto.ShowResult();
}

void LabelBins(TAxis* x) {

    for(Int_t i=1; i<=x->GetNbins(); ++i) {
        const auto a = x->GetBinLowEdge(i);
        x->SetBinLabel(i,to_string(int(a)).c_str());
    }
}

std::shared_ptr<MCReconstructCheck::PositionMap> MCReconstructCheck::histgroup::makePosMap(const HistogramFactory& f,
                                                                                           MCReconstructCheck::histgroup::detectortype d,
                                                                                           const string& name, const string title)
{
    switch (d) {
    case detectortype::All:
    case detectortype::CB:
        return std::make_shared<PositionMapCB>(f,name,title);
    case detectortype::TAPS:
        return std::make_shared<PositionMapTAPS>(f,name,title);
    }
    return nullptr;
}

MCReconstructCheck::histgroup::histgroup(const HistogramFactory& parent, const string& prefix, detectortype d): Prefix(prefix)
{
    HistogramFactory HistFac(prefix, parent, prefix);

    const BinSettings bins_energy(1600);
    const BinSettings bins_vetoEnergy(100,0,10);
    const BinSettings bins_clustersize(19,1,20);
    const BinSettings bins_theta = d == detectortype::TAPS ? BinSettings(50,0,25) : BinSettings(360,0,180);
    const BinSettings bins_anglediff(100,-15,15);

    nPerEvent     = HistFac.makeTH1D("Candidates/Event", "Candidates/Event","",BinSettings(10),"candEvent");
    LabelBins(nPerEvent->GetXaxis());
    nPerEvent->SetFillColor(kGray);

    nPerEventPerE = HistFac.makeTH2D("Candidates/Event/Energy","MC True Energy [MeV]","Candidates/Event",BinSettings(1000),BinSettings(10),"candEventEnergy");
    LabelBins(nPerEventPerE->GetYaxis());

    splitPerEvent = HistFac.makeTH1D("Split-flagged Clusters/Event", "# split clusters/Event","",BinSettings(20),"splits");
    LabelBins(splitPerEvent->GetXaxis());
    splitPerEvent->SetFillColor(kGray);

    splitPos = makePosMap(HistFac,d, "splitpos","Pos Mult > 1");
    splitFlagPos = makePosMap(HistFac,d, "splitflagpos","Pos of split-flagged clusters");
    touchesholeFlagPos = makePosMap(HistFac,d, "touchesholeflagpos","Pos of toucheshole-flagged clusters");
    posCharged = makePosMap(HistFac,d, "chargedpos","Pos of charged cands");

    mult2_split_angles.resize(3);
    mult2_split_stack = HistFac.make<hstack>("Opening Angle Mult==2","mult2_split_stack",true);
    for(size_t i=0;i<mult2_split_angles.size();++i) {
        mult2_split_angles[i] = HistFac.makeTH1D("split="+to_string(i),"#alpha [#circ]","",BinSettings(180,0,90),"mult2_split_"+to_string(i));
        mult2_split_angles[i]->SetLineColor(plot::histstyle::color_t::GetDark(i));
        mult2_split_angles[i]->SetLineWidth(2);
        *mult2_split_stack << mult2_split_angles[i];
    }

    cluserSize = HistFac.makeTH2D("Cluster Size","E [MeV]","Elements",bins_energy, bins_clustersize,"clustersize");
    LabelBins(cluserSize->GetYaxis());

    cluserSize_true = HistFac.makeTH2D("Cluster Size","E_{True} [MeV]","Elements",bins_energy, bins_clustersize,"clustersize_true");
    LabelBins(cluserSize_true->GetYaxis());

    dEE = HistFac.makeTH2D("dEE", "E [MeV]","VetoEnergy [MeV]",bins_energy, bins_vetoEnergy,"dEE");
    dEE_true = HistFac.makeTH2D("dEE (true E)", "E_{True} [MeV]","VetoEnergy [MeV]",bins_energy, bins_vetoEnergy,"dEE_true");

    nCharged        = HistFac.makeTH1D("N Charged (VetoEnergy > 0)", "# charged candidates", "", BinSettings(10),"ncharged");
    LabelBins(nCharged->GetXaxis());
    nCharged->SetFillColor(kGray);

    unmatched_veto = HistFac.makeTH1D("Unmatched Veto Clusters","# unmatched veto clusters","",BinSettings(6),"unmatched_veto");
    LabelBins(unmatched_veto->GetXaxis());
    unmatched_veto->SetFillColor(kGray);

    veto_cand_phi_diff = HistFac.makeTH1D("Angle unmatched Veto - Cand","# unmatched veto clusters","",BinSettings(6),"veto_cand_phi_diff");

    energy_rec_true = HistFac.makeTH2D("Energy","E_{True} [MeV]","E_{Rec} [MeV]",bins_energy,bins_energy,"energy_rec_true");

    theta_diff_truetheta  = HistFac.makeTH2D("Theta Difference","#theta_{True} [#circ]","#theta_{Rec} - #theta_{True} [#circ]",
                                             bins_theta,bins_anglediff,"theta_diff_truetheta");
    phi_diff_truetheta  = HistFac.makeTH2D("Phi Difference","#theta_{True} [#circ]","#phi_{Rec} - #phi_{True} [#circ]",
                                           bins_theta,bins_anglediff,"phi_diff_truetheta");
    openingangle_diff_truetheta  = HistFac.makeTH2D("Opening Angle","#theta_{True} [#circ]","Opening Angle Rec/True [#circ]",
                                                    bins_theta,bins_anglediff,"openingangle_diff_truetheta");

    theta_diff_trueEk  = HistFac.makeTH2D("Theta Difference","E_{True} [MeV]","#theta_{Rec} - #theta_{True} [#circ]",
                                             bins_energy,bins_anglediff,"theta_diff_trueEk");
    phi_diff_trueEk  = HistFac.makeTH2D("Phi Difference","E_{True} [MeV]","#phi_{Rec} - #phi_{True} [#circ]",
                                        bins_energy,bins_anglediff,"phi_diff_trueEk");
    openingangle_diff_trueEk  = HistFac.makeTH2D("Opening Angle","E_{True} [MeV]","Opening Angle Rec/True [#circ]",
                                                 bins_energy,bins_anglediff,"openingangle_diff_trueEk");

    energy_recov  = makePosMap(HistFac,d,"energy_recov","Energy Recovery Average");
    energy_recov->maphist->SetStats(false);

    mc_mult1_positions  = makePosMap(HistFac,d,"mc_mult1_positions","MC position Mult1");
    mc_mult1_positions->maphist->SetStats(false);

    rec_mult1_positions  = makePosMap(HistFac,d,"rec_mult1_positions","Rec position Mult1");
    rec_mult1_positions->maphist->SetStats(false);


    input_positions = makePosMap(HistFac,d,"input","MC True Positions");

    mult1_chargedPos= makePosMap(HistFac,d,"mult1_chargedPos","Mult==1 Charged Pos");


}

void MCReconstructCheck::histgroup::ShowResult() const
{

    canvas c(Prefix);

    c << drawoption("colz") << nPerEvent << nPerEventPerE << splitPerEvent
      << mc_mult1_positions << rec_mult1_positions
      << splitFlagPos << splitPos << touchesholeFlagPos
      << cluserSize << cluserSize_true << dEE << dEE_true << nCharged << posCharged << unmatched_veto
      << drawoption("nostack") << padoption::Legend << mult2_split_stack
      << drawoption("colz") << energy_recov
      << padoption::LogZ << energy_rec_true
      << theta_diff_truetheta << phi_diff_truetheta << openingangle_diff_truetheta
      << theta_diff_trueEk << phi_diff_trueEk << openingangle_diff_trueEk
      << input_positions << mult1_chargedPos
      << endc;
}

void Norm(TH1* hist) {
    hist->Scale(1.0/hist->GetEntries());
}

void MCReconstructCheck::histgroup::Finish()
{
    Norm(nPerEvent);
    Norm(splitPerEvent);
    Norm(nCharged);
    Norm(unmatched_veto);
    energy_recov->maphist->Divide(mc_mult1_positions->maphist);
    mult1_chargedPos->maphist->Divide(mc_mult1_positions->maphist);
}

double angle(const TCandidate& c1, const TCandidate& c2) {
    return vec3::RThetaPhi(1, c1.Theta, c1.Phi).Angle(vec3::RThetaPhi(1, c2.Theta, c2.Phi));
}

TCandidatePtrList CandidatesByDetector(const Detector_t::Any_t& detector, const TCandidateList& candidates) {
    TCandidatePtrList cands;
    for(const auto& c : candidates.get_iter()) {
        if(c->Detector & detector) {
            cands.emplace_back(c);
        }
    }
    return cands;
}


void MCReconstructCheck::histgroup::Fill(const TParticlePtr& mctrue, const TCandidateList& cand, const TClusterList& all_clusters)
{
    const auto mc_phi = std_ext::radian_to_degree(mctrue->Phi());
    const auto mc_theta = mctrue->Theta();
    const auto mc_energy = mctrue->Ek();


    input_positions->Fill(mc_theta, mc_phi);


    nPerEvent->Fill(cand.size());
    nPerEventPerE->Fill(mctrue->Ek(), cand.size());

    unsigned nsplit(0);
    unsigned n(0);
    unsigned ncharged(0);

    for(const TCandidate& c : cand) {
        if(c.VetoEnergy > 0.0) {
            posCharged->Fill(mc_theta,mc_phi);
        }

            cluserSize->Fill(c.CaloEnergy, c.ClusterSize);
            cluserSize_true->Fill(mc_energy, c.ClusterSize);
            dEE->Fill(c.CaloEnergy, c.VetoEnergy);
            dEE_true->Fill(mc_energy, c.VetoEnergy);

            if(c.VetoEnergy > 0.0)
                ncharged++;
            ++n;

        for(const TCluster& cl : c.Clusters) {
            if(cl.HasFlag(TCluster::Flags_t::Split)) {
                ++nsplit;
                splitFlagPos->Fill(mc_theta, mc_phi);
            }
            if(cl.HasFlag(TCluster::Flags_t::TouchesHoleCentral)) {
                touchesholeFlagPos->Fill(mc_theta, mc_phi);
            }
        }
    }

    unsigned nunmatched_veto(0);
    for(const TCluster& ic : all_clusters) {
        if(ic.DetectorType & Detector_t::Any_t::Veto) {
            nunmatched_veto++;
        }
    }
    unmatched_veto->Fill(nunmatched_veto);

    if(n>0)
        nCharged->Fill(ncharged);

    splitPerEvent->Fill(nsplit);

    if(cand.size() == 2) {
        if(nsplit<3) {
            mult2_split_angles[nsplit]->Fill(std_ext::radian_to_degree(angle(cand.at(0), cand.at(1))));
        }
    }

    if(cand.size()==1) {
        const auto& c = cand.at(0);
        const auto rec = c.CaloEnergy / mc_energy;
        energy_recov->Fill(mc_theta, mc_phi, rec);
        mc_mult1_positions->Fill(mc_theta,mc_phi);
        rec_mult1_positions->Fill(c.Theta,std_ext::radian_to_degree(c.Phi));
        energy_rec_true->Fill(mc_energy,c.CaloEnergy);

        theta_diff_truetheta->Fill(std_ext::radian_to_degree(mc_theta), std_ext::radian_to_degree(c.Theta - mc_theta));
        phi_diff_truetheta->Fill(std_ext::radian_to_degree(mc_theta), std_ext::radian_to_degree(c.Phi) - mc_phi);
        openingangle_diff_truetheta->Fill(std_ext::radian_to_degree(mc_theta),std_ext::radian_to_degree(mctrue->Angle(c)));

        theta_diff_trueEk->Fill(mc_energy, std_ext::radian_to_degree(c.Theta - mc_theta));
        phi_diff_trueEk->Fill(mc_energy, std_ext::radian_to_degree(c.Phi) - mc_phi);
        openingangle_diff_trueEk->Fill(mc_energy,std_ext::radian_to_degree(mctrue->Angle(c)));

        if(c.VetoEnergy>0.0) {
            mult1_chargedPos->Fill(mc_theta,mc_phi);
        }
    } else if( cand.size() > 1) {
        splitPos->Fill(mc_theta,mc_phi);
    }

}



void MCReconstructCheck::PositionMapCB::Fill(const double theta, const double phi, const double v)
{
    maphist->Fill(cos(theta),phi,v);
}

void MCReconstructCheck::PositionMap::Draw(const string&) const
{
    maphist->Draw("colz");
}

MCReconstructCheck::PositionMapCB::PositionMapCB(const HistogramFactory &f, const string &name, const string &title)
{
    const BinSettings costheta(360,-1,1);
    const BinSettings phi(360,-180,180);
    maphist = f.makeTH2D(title,"cos(#theta_{True})","#phi [#circ]",costheta,phi,name);
}

MCReconstructCheck::PositionMapTAPS::PositionMapTAPS(const HistogramFactory &f, const string &name, const string& title) :
    taps(ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>()),
    taps_dist(taps->GetZPosition())
{
    const auto radius = taps->GetRadius();

    const BinSettings bins(2*radius/0.5,-radius,radius);

    maphist = f.makeTH2D(title,"x [cm]","y [cm]",bins, bins, name);
}

void MCReconstructCheck::PositionMapTAPS::Fill(const double theta, const double phi, const double v)
{
    constexpr const auto max_theta = std_ext::degree_to_radian(25.0);

    if(theta < max_theta) {
        const auto r = tan(theta) * taps_dist; //cm
        const auto x = r * cos(phi);
        const auto y = r * sin(phi);
        maphist->Fill(x,y,v);
    }
}

void MCReconstructCheck::PositionMapTAPS::Draw(const string&) const
{
    maphist->Draw("colz");
    auto taps = new TH2TAPS();
    taps->Draw("same text");
}



MCReconstructCheck::TAPSVetoMatch::TAPSVetoMatch(const HistogramFactory& parent)
{
    HistogramFactory f("TapsVetoMatch", parent);
    auto tapsveto = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPSVeto);
    vetoElement_dist = f.makeTH2D("Veto TAPS Cluster dist","Veto Element","Dist [cm]",
                                  BinSettings(tapsveto->GetNChannels()),BinSettings(100,0,20),"cluster_dist");
}

void MCReconstructCheck::TAPSVetoMatch::ShowResult()
{
    canvas c("TAPS Veto");
    c << vetoElement_dist << endc;
}

void MCReconstructCheck::TAPSVetoMatch::Fill(const TCandidateList& cands, const TClusterList& all_clusters)
{
    using namespace ant::std_ext;

    auto clusterLoop = [this] (const TClusterPtr& vCluster, const TCandidateList& cands) {
        if(vCluster && vCluster->DetectorType == Detector_t::Type_t::TAPSVeto) {
            for(const TCandidate& cCand : cands) {
                const auto cCluster = cCand.FindCaloCluster();
                if(cCluster && cCluster->DetectorType == Detector_t::Type_t::TAPS) {
                    const auto dx = vCluster->Position.x - cCluster->Position.x;
                    const auto dy = vCluster->Position.y - cCluster->Position.y;
                    const auto d  = sqrt(sqr(dx)+sqr(dy));
                    this->vetoElement_dist->Fill(vCluster->CentralElement, d);
                }
            }
        }
    };


    for(const TCandidate& vCand : cands) {
        const auto vCluster = vCand.FindVetoCluster();
        clusterLoop(vCluster, cands);
    }

    for(auto iCluster=all_clusters.begin(); iCluster != all_clusters.end(); ++iCluster) {
        clusterLoop(iCluster.get_ptr(), cands);
    }

}

AUTO_REGISTER_PHYSICS(MCReconstructCheck)
