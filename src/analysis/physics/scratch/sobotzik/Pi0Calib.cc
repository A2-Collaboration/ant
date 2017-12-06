#include "Pi0Calib.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/Logger.h"
#include "TH1D.h"
#include "TTree.h"
#include "vector"
#include "TStyle.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

scratch_sobotzik_Pi0Calib::scratch_sobotzik_Pi0Calib(const string& name, OptionsPtr opts)
    : Physics(name, opts),
      CaloEnergy_Window(opts->Get<interval<double>>("CaloEnergyWindow", {-std_ext::inf, std_ext::inf})),
      promptrandom(ExpConfig::Setup::Get())
{


    const BinSettings bins_IM   (400, 0, 1100); // MeV
    const BinSettings bins_angle(180, 0,    180); // degrees
    const BinSettings bins_timing(200,-100,100); // ns
    const BinSettings bins_energy(200, 0, 1000);

    h_IM_All   = HistFac.makeTH1D("IM: All",  "IM / MeV","",bins_IM,"IM_All");

    h_CB_E_True_Opening_Angle = HistFac.makeTH2D("IM: CB true Opening Angle & Rec. Energy", "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_True_Angle");
    h_CB_Angle_True_E_Angle = HistFac.makeTH2D("IM: CB true Energy & Rec. Opening Angle", "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_True_E");

    h_IM_CB_ClustersizeOAngle     = HistFac.makeTH2D("IM: CB Clustersize vs. OAngle","Opening Angle [^{#circ}]","Clustersize",BinSettings(180,0,180),BinSettings(20,0,20),"IM_CB_ClusterOAngle");

    h_IM_CB_all             = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_All");

    h_Meson_Energy_interval =HistFac.makeTH3D("MC-Meson-Symmetric-Photons","IM / MeV", "E_{#gamma} [MeV]", "Meson Energy [MeV]",bins_IM,BinSettings(32,0,800),BinSettings(158,0,1580),"Meson_Energy_Interval");
    h_Meson_Energy_interval_30_Degree_Cut =HistFac.makeTH3D("MC-Meson-Symmetric-Photons","IM / MeV", "E_{#gamma} [MeV]", "Meson Energy [MeV]",bins_IM,BinSettings(32,0,800),BinSettings(158,0,1580),"Meson_Energy_Interval_30_Degree_Cut");

    h_IM_CB_Uncharged_No_Cut             = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Uncharged");
    h_IM_CB_interval        = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Interval");
    h_IM_CB_interval_Uncharged_No_Cut        = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Interval_No_Cut");
    h_IM_CB_interval_Uncharged_30_Degree_Cut        = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Interval_30_Degree_Cut");

//    h_IM_CB_One_high_Photon = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_One_high_Photon");

    h_IM_CB_Uncharged_30_Degree_Cut    = HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Uncharged_30_Degree_Cut");


    h_IM_CB_Angle_Energy    = HistFac.makeTH2D("IM: Angle",   "Angle / Degrees","E_{#gamma} [MeV]",bins_angle,BinSettings(32,0,800),"IM_CB_Angle");

    h_IM_CB_Min_Opening_Angle =HistFac.makeTH2D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Min_Opening_Angle");
    h_IM_CB_Rec_vs_Gen_Opening_Angle = HistFac.makeTH3D("Rec. vs. Gen. Opening Angle","Reconstructed Opening Angle / Degree", "Generated Opening Angle / Degree","Energ of the Photons E[MeV]",BinSettings(180,0,180),BinSettings(180,0,180),BinSettings(32,0,800),"IM_CB_Rec_vs_Gen_Opening_Angle");
    h_IM_CB_Rec_vs_Gen_Opening_Angle_Deviation = HistFac.makeTH2D("IM: Deviation between Gen. and Rec. Opening Angle",   "Angle / Degrees","E_{#gamma} [MeV]",BinSettings(200,-10,10),BinSettings(32,0,800),"IM_CB_Rec_vs_Gen_Opening_Angle_Deviation");

    //    h_IM_CB_Rec_vs_Gen_Energie = HistFac.makeTH2D ("Rec. vs. Gen. Energy", "Reconstructed Energy [MeV]" , "Generated Energy [MeV]", BinSettings(1000,0,1000), BinSettings(1000,0,1000),"IM_CB_Rec_vs_Gen_Energy" );
    //    h_IM_CB_Rec_Gen_Energie_Deviation= HistFac.makeTH2D ("E(rec) - E(gen)", "Deviation of the energies [MeV]" , "Energy of the detected Photons E_{#gamma} [MeV]", BinSettings(80,-40,40), BinSettings(32,0,800),"IM_CB_Deviation_Gen_Rec_Energy" );


    //    h_IM_CB_Theta_Phi_Energy= histFac.makeTH3D("IM:CB","Polar angle Theta / Degree","Azimut angle Phi / Degree","Energy of the Photons E_{#gamma} [MeV]", bins_angle,BinSettings(360,-180,180) ,BinSettings(32,0,800),"IM_CB_Theta_Phi");
    h_IM_CB_interval_Theta_Phi_Energy= HistFac.makeTH3D("IM:CB","Polar angle Theta / Degree","Azimut angle Phi / Degree","Energy of the Photons E_{#gamma} [MeV]", bins_angle,BinSettings(360,-180,180) ,BinSettings(32,0,800),"IM_CB_Interval_Theta_Phi");

    //    h_IM_CB_ZVertex         = histFac.makeTH3D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]","Z-Vertex [cm]",bins_IM,BinSettings(32,0,800),BinSettings(10,-5,5),"IM_CB_ZVertex");
    //    h_IM_CB_ZVertex_interval         = histFac.makeTH3D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]","Z-Vertex [cm]",bins_IM,BinSettings(32,0,800),BinSettings(10,-5,5),"IM_CB_ZVertex_interval");
    h_IM_CB_ZVertex_interval_30_Degree_Cut         = HistFac.makeTH3D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]","Z-Vertex [cm]",bins_IM,BinSettings(32,0,800),BinSettings(10,-5,5),"IM_CB_ZVertex_interval_30_Degree_Cut");
    h_IM_CB_ZVertex         = HistFac.makeTH3D("IM: CB",   "IM / MeV","E_{#gamma} [MeV]","Z-Vertex [cm]",bins_IM,BinSettings(32,0,800),BinSettings(10,-5,5),"IM_CB_ZVertex");

    h_IM_CB_AngleDeviation_Energy   = HistFac.makeTH2D("IM: Angle Deviation between Gen. and rec. Photons",   "Angle / Degrees","E_{#gamma} [MeV]",BinSettings(20,0,20),BinSettings(32,0,800),"IM_CB_AngleDeviation");
//    h_IM_CB_AngleDeviation_Photon_Meson_Energy = HistFac.makeTH3D("IM: CB",   "IM / MeV", "Deviation of the opening angle in Degree","Meson Energy [MeV]",bins_IM,BinSettings(20,0,20),BinSettings(158,0,1580),"IM_CB_AngleDeviation_Meson");

    h_IM_CB_corr    = HistFac.makeTH1D("IM: CB corr",   "IM / MeV","",bins_IM,"IM_CB_corr");
    h_IM_TAPS  = HistFac.makeTH1D("IM: TAPS", "IM / MeV","",bins_IM,"IM_TAPS");

    h_Angle_CB   = HistFac.makeTH1D("Angle: CB",   "angle [#circ]","",bins_angle,"Angle_CB");
    h_Angle_TAPS = HistFac.makeTH1D("Angle: TAPS", "angle [#circ]","",bins_angle,"Angle_TAPS");

    h_IM_True_Opening_Angle =HistFac.makeTH2D("True Opening Angle", "Opening angle #alpha / Deg","E_{#gamma} [MeV]",BinSettings(360,0,180),BinSettings(32,0,800),"IM_True_OpeningAngle" );
    h_IM_Rec_Opening_Angle =HistFac.makeTH2D("Rec Opening Angle", "Opening angle #alpha / Deg","E_{#gamma} [MeV]",BinSettings(360,0,180),BinSettings(32,0,800),"IM_Rec_OpeningAngle" );



    h_ClusterHitTiming_CB   = HistFac.makeTH2D("ClusterHitTiming: CB",   "Energy","t / ns",bins_energy,bins_timing,"ClusterHitTiming_CB");
    h_ClusterHitTiming_TAPS = HistFac.makeTH2D("ClusterHitTiming: TAPS", "Energy","t / ns",bins_energy,bins_timing,"ClusterHitTiming_TAPS");

    h_IM_CB_ClusterSize3 = HistFac.makeTH2D("IM Clustersize > 3", "IM / MeV","E_{#gamma} [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_ClusterSize3");

    h_IM_CB_InvOAngletrue = HistFac.makeTH2D("True OpeningAngle vs. Inv. Mass","IM /MeV","True Opening Angle [^{#circ}]",bins_IM,BinSettings(180,0,180),"IM_CB_InvOAngletrue");
    h_IM_CB_InvOAnglerec = HistFac.makeTH2D("Rec OpeningAngle vs. Inv. Mass","IM /MeV","Rec Opening Angle [^{#circ}]",bins_IM,BinSettings(180,0,180),"IM_CB_InvOAnglerec");
    h_IM_CB_NClusterEnergy= HistFac.makeTH2D("Number of Clusters vs. Energy","E_{#gamma} [MeV]","Number of Clusters",bins_energy,BinSettings(50,0,50),"IM_CB_NCluster");

//    h_CB_Theta_Diff = HistFac.makeTH3D("#Theta_{true} - #Theta_{rec} vs. #Theta_{rec} for different energies","#Theta_{rec} [#circ]","#Theta_{true} - #Theta_{rec}","E_{#gamma} [MeV]",BinSettings(180,0,180),BinSettings(20,-10,10),BinSettings(32,0,800),"CB_Theta_Diff");
//    h_CB_Theta_Diff = HistFac.makeTH3D("#Phi_{true} - #Phi_{rec} vs. #Phi_{rec} for different energies","#Phi_{rec} [#circ]","#Phi_{true} - #Phi_{rec}","E_{#gamma} [MeV]",BinSettings(180,0,180),BinSettings(20,-10,10),BinSettings(32,0,800),"CB_Phi_Diff");








    for(int i=0;i<8;++i) {
        const string name_Symmetric = std_ext::formatter() << "CB " << i * 100 <<" MeV to "<<(i+1) * 100<<" MeV Clustersize > 3";

        h_cbs_ClusterSize3.push_back(HistFac.make<TH2CB>(name_Symmetric.c_str(),name_Symmetric.c_str()));

    }
    for( int i = 0; i< 8; ++i)
    {
        const string name_All_Photons = std_ext::formatter() << "CB " << i * 100 <<" MeV to "<<(i+1) * 100<<" MeV Clustersize > 0";
        h_cbs_ClusterSize0.push_back(HistFac.make<TH2CB>(name_All_Photons.c_str(),name_All_Photons.c_str()));
    }


    t.CreateBranches(HistFac.makeTTree("cluster_sym"));

}




TParticleTree_t getFirst(const ParticleTypeDatabase::Type& t, const TParticleTree_t& tree) {
    auto node = tree->Get();
    if(node->Type() == t) {
        return tree;
    } else {
        for(const auto& d : tree->Daughters()){
            auto r = getFirst(t,d);
            if(r)
                return r;
        }
    }
    return nullptr;
}

void scratch_sobotzik_Pi0Calib::ProcessEvent(const TEvent& event, manager_t&)
{


    auto ptree = event.MCTrue().ParticleTree;
    TParticleTree_t true_pi0_tree = nullptr;
    if(ptree) {
        //        auto typetree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g);
        //        if(!ptree->IsEqual(typetree, utils::ParticleTools::MatchByParticleName))
        //            return;
        true_pi0_tree = getFirst(ParticleTypeDatabase::Pi0, ptree);
        if(!true_pi0_tree)
            return;
    }

    TCandidatePtrList c_CB;
    TCandidatePtrList c_TAPS;
    for(auto& c : event.Reconstructed().Candidates.get_iter()) {
        if(c->Detector & Detector_t::Type_t::CB)
            c_CB.emplace_back(c);
        else if(c->Detector & Detector_t::Type_t::TAPS)
            c_TAPS.emplace_back(c);
    }
    auto zVertex = event.MCTrue().Target.Vertex.z;

    TParticlePtr true_pi0 = nullptr;
    vector<TParticlePtr> true_gamma;

    if(true_pi0_tree) {
        true_pi0 = true_pi0_tree->Get();
        if(true_pi0_tree->Daughters().size() == 2) {
            true_gamma.push_back(true_pi0_tree->Daughters().front()->Get());
            true_gamma.push_back(true_pi0_tree->Daughters().back()->Get());
        }
    }



    if(c_CB.size() != 2)
        return;



    auto sum_as_photons = [this] (const TCandidatePtrList& cands) {
        LorentzVec sum;
        for(auto& cand : cands) {
            if (CaloEnergy_Window.Contains(cand-> CaloEnergy))
            {
                sum += TParticle(ParticleTypeDatabase::Photon, cand);
            }
        }
        return sum;
    };

    auto sum_as_corr_photons = [] (const TCandidatePtrList& cands) {
        // copied from Sergey's kinfitter header
        // to be applied as Ecorr = Ecl*(1+Fcor)
        auto EclCorCB = [] (Double_t Ecl) {
            // CB energy correction
            // correction for 12 MeV cluster threshold in the CB.
            Double_t p[5] = {1.52950e-02, 5.92991e-03, 4.57857e-01, 8.98180e-03,
                             7.75002e-03}; // photon smeared smcal11
            Double_t Fcor = p[0] / pow(Ecl + p[1], p[2]) + p[3] + p[4] * Ecl;
            return Fcor;
        };

        LorentzVec sum;
        for(auto& cand : cands) {
            TParticle p(ParticleTypeDatabase::Photon, cand);
            const double corr = EclCorCB(p.Ek()/1000.0);
            p *= 1+corr;
            sum += p;
        }
        return sum;

    };





    const auto min_angle = [] (const TCandidatePtrList& cands) {

        double angle = std_ext::inf;

        for(auto c = utils::makeCombination(cands,2); !c.done(); ++c) {
            vec3 c0(*c.at(0));
            auto& c1 = *c.at(1);
            angle = min(angle, std_ext::radian_to_degree(c0.Angle(c1)));
        };
        return angle;
    };

    const auto fill_timing = [] (const TCandidatePtrList& cands, TH2D* h) {
        for(auto& cand : cands) {
            auto cl = cand->FindCaloCluster();

            for(auto& hit : cl->Hits)
                h->Fill(hit.Energy, hit.Time);
            //            h->Fill(cl->Energy, cl->Time);
        }
    };




    std::array<double,2> min_angle_rg;
    std::array<int,2> j;
    std::array<double,2> true_gamma_energy;
    int iter = 0;
    int clen = c_CB.size();
    double  rec_opening_angle  = 0;
    double  true_opening_angle = 0;
    auto true_theta1 =0.0;
    auto true_theta2 =0.0;
    auto true_phi1 =0.0;
    auto true_phi2 =0.0;
    if(true_pi0_tree)
    {
        for(const auto& gamma : true_gamma)
        {
            double min_angle= std_ext::inf;
            for(int n=0; n<clen; n++)
            {
                if(static_cast<vec3>(*c_CB.at(n)).Angle(gamma->p) < min_angle)
                {
                    min_angle = static_cast<vec3>(*c_CB.at(n)).Angle(gamma->p);
                    j[iter]=n;
                }
            }
            min_angle_rg[iter] = min_angle;
            true_gamma_energy[iter] = gamma->Ek();
            iter++;
        };


        std::array<vec3,2> true_cand_array;
        int accumulator = 0;
        for(const auto& gamma : true_gamma)
        {
            true_cand_array[accumulator] = gamma->p;
            accumulator++;
        }

        //Calculation of True and reconstructed opening angle
        //opening_angle between the candidates
        rec_opening_angle  = static_cast<vec3>(*c_CB.at(0)).Angle(*c_CB.at(1));

        true_opening_angle = true_cand_array[0].Angle(true_cand_array[1]);

        true_theta1 = true_cand_array[0].Theta();
        true_theta2 = true_cand_array[1].Theta();
        true_phi1   = true_cand_array[0].Phi();
        true_phi2   = true_cand_array[1].Phi();

    }

    double angleedge = 30;
    const auto& sum_CB = sum_as_photons(c_CB);
    const auto& sum_TAPS = sum_as_photons(c_TAPS);
    const auto& angle_CB = min_angle(c_CB);

    triggersimu.ProcessEvent(event);
    auto w = event.Reconstructed().TaggerHits.size() == 0 ? 1.0 : 0.0;
    for(const TTaggerHit& TagH : event.Reconstructed().TaggerHits) {

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(TagH));

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        w += promptrandom.FillWeight();
    }

    h_IM_All->Fill((sum_CB+sum_TAPS).M(), w);

    // only symmetric Photons
    const auto binwidth = h_IM_CB_interval->GetYaxis()->GetBinWidth(1);
    //    const auto bin2 = h_IM_CB_interval->GetYaxis()->FindBin(c_CB.at(1)->CaloEnergy);
    auto bindiff= c_CB.at(0)->CaloEnergy - c_CB.at(1)->CaloEnergy;
    if(bindiff<0){
        bindiff *= -1;
    }

    if(bindiff <= binwidth) {
        if(sum_CB.M()>1.0) {


            t.E1 = c_CB.at(0)->CaloEnergy;
            t.E2 = c_CB.at(1)->CaloEnergy;
            t.M  = sum_CB.M();
            t.Theta1_rec = c_CB.at(0)->Theta;
            t.Theta2_rec = c_CB.at(1)->Theta;
            t.Phi1_rec   = c_CB.at(0)->Phi;
            t.Phi2_rec   = c_CB.at(1)->Phi;

            t.ClusterSize1 = c_CB.at(0)->FindCaloCluster()->Hits.size();
            t.ClusterSize2 = c_CB.at(1)->FindCaloCluster()->Hits.size();
            t.OpeningAngle = rec_opening_angle;
            t.ClusterNumber1 = c_CB.at(0)->FindCaloCluster()->CentralElement;
            t.ClusterNumber2 = c_CB.at(1)->FindCaloCluster()->CentralElement;
            t.w = promptrandom.FillWeight();



            if(true_pi0){
                t.ZVertex = zVertex;
                t.true_E1 = true_gamma_energy[0];
                t.true_E2 = true_gamma_energy[1];
                t.true_openingangle = true_opening_angle;
                t.true_m = sqrt(2 * true_gamma_energy[0] * true_gamma_energy[1] * (1 - cos(true_opening_angle)));
                t.Theta1_true=true_theta1;
                t.Theta2_true=true_theta2;
                t.Phi1_true  =true_phi1;
                t.Phi2_true  =true_phi2;
                t.Cand1_Angle_rectrue=min_angle_rg[0];
                t.Cand2_Angle_rectrue=min_angle_rg[1];

            }

            LOG(INFO)<<true_theta1<<"Rec:"<< c_CB.at(0)->Theta <<endl;
            t.Tree->Fill();



            h_IM_CB_ClustersizeOAngle->Fill(std_ext::radian_to_degree(rec_opening_angle),c_CB.at(0)->FindCaloCluster()->Hits.size());
            h_IM_CB_ClustersizeOAngle->Fill(std_ext::radian_to_degree(rec_opening_angle),c_CB.at(1)->FindCaloCluster()->Hits.size());


            //Cut on pi0
            if(sum_CB.M() > 70.0 && sum_CB.M() < 220.0)
            {
                const auto cluster1 = c_CB.at(0)->FindCaloCluster();
                const auto cluster2 = c_CB.at(1)->FindCaloCluster();

                if(cluster1 && cluster2)
                {
                    if(cluster1->Hits.size() > 0 && cluster2->Hits.size() > 0)
                    {

                        int j1  = c_CB.at(0)->CaloEnergy / 100.0;
                        int j2  = c_CB.at(1)->CaloEnergy / 100.0;
                        if (j1 < 8 ){

                            h_cbs_ClusterSize0.at(j1)->FillElement(cluster1->CentralElement,1);
                        }

                        if   ( j2 < 8 ){
                            h_cbs_ClusterSize0.at(j2)->FillElement(cluster2->CentralElement,1);
                        }

                    }
                }


                if(cluster1 && cluster2)
                {
                    if(cluster1->Hits.size() > 3 && cluster2->Hits.size() > 3)
                    {


                        int j1  = c_CB.at(0)->CaloEnergy / 100.0;
                        int j2  = c_CB.at(1)->CaloEnergy / 100.0;

                        if( j1 < 8 && j2 < 8){
                            h_cbs_ClusterSize3.at(j1)->FillElement(cluster1->CentralElement,1);
                            h_cbs_ClusterSize3.at(j2)->FillElement(cluster2->CentralElement,1);

                            h_IM_CB_ClusterSize3->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
                            h_IM_CB_ClusterSize3->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy);
                        }
                    }
                }

            }


            h_IM_True_Opening_Angle->Fill(std_ext::radian_to_degree(true_opening_angle),true_gamma_energy[0]);
            h_IM_True_Opening_Angle->Fill(std_ext::radian_to_degree(true_opening_angle),true_gamma_energy[1]);
            h_IM_Rec_Opening_Angle->Fill(std_ext::radian_to_degree(rec_opening_angle),true_gamma_energy[0]);
            h_IM_Rec_Opening_Angle->Fill(std_ext::radian_to_degree(rec_opening_angle),true_gamma_energy[1]);




            h_IM_CB_interval->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,w);
            h_IM_CB_interval->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,w);

            //            h_IM_CB_ZVertex_interval->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,zVertex);
            //            h_IM_CB_ZVertex_interval->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,zVertex);

            if((c_CB.at(0)->VetoEnergy == 0 )&&(c_CB.at(1)->VetoEnergy == 0))
            {
                if(true_pi0) {
                    h_Meson_Energy_interval->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,true_pi0->Ek());
                    h_Meson_Energy_interval->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,true_pi0->Ek());
                }

                h_IM_CB_interval_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,w);
                h_IM_CB_interval_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,w);


                if(sum_CB.M() > 70.0 && sum_CB.M() < 220.0)
                {
                    h_IM_CB_interval_Theta_Phi_Energy->Fill(c_CB.at(0)->Theta / (2 * 3.141) *360,c_CB.at(0)->Phi / (2 * 3.141) *360, c_CB.at(0)->CaloEnergy);
                    h_IM_CB_interval_Theta_Phi_Energy->Fill(c_CB.at(1)->Theta / (2 * 3.141) *360,c_CB.at(1)->Phi / (2 * 3.141) *360 ,c_CB.at(1)->CaloEnergy);
                }

                if(true_pi0)
                {
                    //                    h_IM_CB_Rec_vs_Gen_Energie->Fill(c_CB.at(j[0])->CaloEnergy,true_gamma_energy[0]);
                    //                    h_IM_CB_Rec_vs_Gen_Energie->Fill(c_CB.at(j[1])->CaloEnergy,true_gamma_energy[1]);

                    //                    h_IM_CB_Rec_Gen_Energie_Deviation->Fill((c_CB.at(j[0])->CaloEnergy-true_gamma_energy[0]),c_CB.at(0)->CaloEnergy);
                    //                    h_IM_CB_Rec_Gen_Energie_Deviation->Fill((c_CB.at(j[1])->CaloEnergy-true_gamma_energy[1]),c_CB.at(1)->CaloEnergy);
                }

                //checking the opening angle between the candidates; only fill if the angle is 30 Degree or higher
                if(rec_opening_angle > std_ext::degree_to_radian(30.0))
                {
                    h_IM_CB_Min_Opening_Angle->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,w);
                    h_IM_CB_Min_Opening_Angle->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,w);

                }

                if(true_pi0)
                {
                    h_IM_CB_Rec_vs_Gen_Opening_Angle->Fill(std_ext::radian_to_degree(rec_opening_angle),std_ext::radian_to_degree(true_opening_angle),c_CB.at(0)->CaloEnergy);
                    h_IM_CB_Rec_vs_Gen_Opening_Angle->Fill(std_ext::radian_to_degree(rec_opening_angle),std_ext::radian_to_degree(true_opening_angle),c_CB.at(1)->CaloEnergy);
                    h_IM_CB_Rec_vs_Gen_Opening_Angle_Deviation->Fill(std_ext::radian_to_degree(rec_opening_angle) - std_ext::radian_to_degree(true_opening_angle),c_CB.at(0)->CaloEnergy);
                    h_IM_CB_Rec_vs_Gen_Opening_Angle_Deviation->Fill(std_ext::radian_to_degree(rec_opening_angle) - std_ext::radian_to_degree(true_opening_angle),c_CB.at(1)->CaloEnergy);

                    h_CB_E_True_Opening_Angle->Fill(sqrt(2 * c_CB.at(0)->CaloEnergy * c_CB.at(1)->CaloEnergy * (1-cos(true_opening_angle))),c_CB.at(0)->CaloEnergy);
                    h_CB_E_True_Opening_Angle->Fill(sqrt(2 * c_CB.at(0)->CaloEnergy * c_CB.at(1)->CaloEnergy * (1-cos(true_opening_angle))),c_CB.at(1)->CaloEnergy);

                    h_CB_Angle_True_E_Angle->Fill(sqrt(2 * true_gamma_energy[0] * true_gamma_energy[1] * (1-cos(rec_opening_angle))),c_CB.at(0)->CaloEnergy);
                    h_CB_Angle_True_E_Angle->Fill(sqrt(2 * true_gamma_energy[1] * true_gamma_energy[1] * (1-cos(rec_opening_angle))),c_CB.at(1)->CaloEnergy);

                }



                if(     (c_CB.at(0)->Theta >(angleedge * 2 * 3.141 /360) &&
                         c_CB.at(0)->Theta <180 - (angleedge * 2 * 3.141 /360))
                        &&
                        (c_CB.at(1)->Theta >(angleedge * 2 * 3.141 /360) &&
                         c_CB.at(1)->Theta <180 - (angleedge * 2 * 3.141 /360)))


                {
                    if(true_pi0) {
                        h_Meson_Energy_interval_30_Degree_Cut->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,true_pi0->Ek());
                        h_Meson_Energy_interval_30_Degree_Cut->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,true_pi0->Ek());


                        h_IM_CB_AngleDeviation_Energy->Fill(std_ext::radian_to_degree(min_angle_rg[0]), c_CB.at(j[0])-> CaloEnergy);
                        h_IM_CB_AngleDeviation_Energy->Fill(std_ext::radian_to_degree(min_angle_rg[1]), c_CB.at(j[1])-> CaloEnergy);

//                        h_IM_CB_AngleDeviation_Photon_Meson_Energy->Fill(std_ext::radian_to_degree(min_angle_rg[0]),c_CB.at(j[0])-> CaloEnergy,true_pi0->Ek());
//                        h_IM_CB_AngleDeviation_Photon_Meson_Energy->Fill(std_ext::radian_to_degree(min_angle_rg[1]),c_CB.at(j[1])-> CaloEnergy,true_pi0->Ek());

                    }

                    h_IM_CB_interval_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(0)->CaloEnergy,w);
                    h_IM_CB_interval_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(1)->CaloEnergy,w);

                    h_IM_CB_ZVertex_interval_30_Degree_Cut->Fill (sum_CB.M(),c_CB.at(0)->CaloEnergy,zVertex,w);
                    h_IM_CB_ZVertex_interval_30_Degree_Cut->Fill (sum_CB.M(),c_CB.at(1)->CaloEnergy,zVertex,w);
                }
            }

        }
    }

    //All Photons allowed

    h_IM_CB_NClusterEnergy->Fill(c_CB.at(0)->CaloEnergy,c_CB.at(0)->FindCaloCluster()->Hits.size());
    h_IM_CB_NClusterEnergy->Fill(c_CB.at(1)->CaloEnergy,c_CB.at(1)->FindCaloCluster()->Hits.size());
    if(sum_CB.M()>1.0)
    {
        t.All_E1_rec = c_CB.at(0)->CaloEnergy;
        t.All_E2_rec = c_CB.at(1)->CaloEnergy;
        t.All_OpeningAngle_rec = rec_opening_angle;
        t.All_Phi1_rec= c_CB.at(0)->Phi;
        t.All_Phi2_rec= c_CB.at(1)->Phi;
        t.All_Theta1_rec= c_CB.at(0)->Theta;
        t.All_Theta2_rec= c_CB.at(1)->Theta;
        t.All_M=sum_CB.M();





        if(true_pi0){
            t.All_E1_true = true_gamma_energy[0];
            t.All_E2_true = true_gamma_energy[1];
            t.All_Phi1_true=true_phi1;
            t.All_Phi2_true=true_phi2;
            t.All_Theta1_true=true_theta1;
            t.All_Theta2_true=true_theta2;
            t.All_OpeningAngle_true= true_opening_angle;
            t.All_Cand1_Angle_rectrue = min_angle_rg[0];
            t.All_Cand2_Angle_rectrue = min_angle_rg[1];
        }



        h_IM_CB_all->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
        h_IM_CB_all->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy);
        h_IM_CB_Angle_Energy->Fill( angle_CB,c_CB.at(0)->CaloEnergy);
        h_IM_CB_Angle_Energy->Fill( angle_CB,c_CB.at(1)->CaloEnergy);
        h_IM_CB_InvOAnglerec->Fill(sum_CB.M(),std_ext::radian_to_degree(true_opening_angle));

        if(true_pi0){
                    h_IM_CB_InvOAngletrue ->Fill(sum_CB.M(),std_ext::radian_to_degree(rec_opening_angle));
        }


        //        h_IM_CB_ZVertex->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,zVertex);
        //        h_IM_CB_ZVertex->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,zVertex);


        if((c_CB.at(0)->VetoEnergy == 0 )&&(c_CB.at(1)->VetoEnergy == 0))
        {
            h_IM_CB_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
            h_IM_CB_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy);
            //            h_IM_CB_Theta_Phi_Energy->Fill(c_CB.at(0)->Theta / (2 * 3.141) *360,c_CB.at(0)->Phi / (2 * 3.141) *360, c_CB.at(0)->CaloEnergy);
            //            h_IM_CB_Theta_Phi_Energy->Fill(c_CB.at(1)->Theta / (2 * 3.141) *360,c_CB.at(1)->Phi / (2 * 3.141) *360 ,c_CB.at(1)->CaloEnergy);

            if(     (c_CB.at(0)->Theta >(angleedge * 2 * 3.141 /360) &&
                     c_CB.at(0)->Theta <180 - (angleedge * 2 * 3.141 /360))
                    &&
                    (c_CB.at(1)->Theta >(angleedge * 2 * 3.141 /360) &&
                     c_CB.at(1)->Theta <180 - (angleedge * 2 * 3.141 /360)))
            {
                h_IM_CB_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(0)->CaloEnergy);
                h_IM_CB_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(1)->CaloEnergy);

//                if(c_CB.at(0)->CaloEnergy > c_CB.at(1)->CaloEnergy)
//                {
//                    h_IM_CB_One_high_Photon->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
//                }
//                else
//                {
//                    h_IM_CB_One_high_Photon->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy);
//                }
            }

            h_IM_CB_ZVertex->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,zVertex,w);
            h_IM_CB_ZVertex->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy,zVertex,w);
        }


    }

    h_IM_CB_corr->Fill(sum_as_corr_photons(c_CB).M());
    h_IM_TAPS->Fill(sum_TAPS.M());

    h_Angle_CB->Fill(min_angle(c_CB));
    h_Angle_TAPS->Fill(min_angle(c_TAPS));


    fill_timing(c_CB, h_ClusterHitTiming_CB);
    fill_timing(c_TAPS, h_ClusterHitTiming_TAPS);
    gStyle->SetOptStat(0);

}

void scratch_sobotzik_Pi0Calib::ShowResult()
{
    gStyle->SetOptStat(0);

    canvas c(GetName());
    //            << h_Angle_CB
    //            << h_Angle_TAPS
    //            << h_IM_All

    c << drawoption("colz")

      << h_IM_CB_all
      << h_IM_CB_interval
      << h_IM_CB_NClusterEnergy
      << h_IM_CB_ClustersizeOAngle
      << h_IM_CB_InvOAngletrue
      << h_IM_CB_InvOAnglerec
      << h_CB_E_True_Opening_Angle
      << h_CB_Angle_True_E_Angle
      << h_IM_CB_interval_Uncharged_No_Cut
      << h_IM_CB_interval_Uncharged_30_Degree_Cut
      << h_IM_CB_Uncharged_No_Cut
      << h_IM_CB_Angle_Energy
         //            << h_IM_CB_Theta_Phi_Energy
      << h_IM_CB_Min_Opening_Angle
      << h_IM_CB_Rec_vs_Gen_Opening_Angle
      << h_IM_CB_Rec_vs_Gen_Opening_Angle_Deviation
         //            << h_IM_CB_Rec_vs_Gen_Energie
         //            << h_IM_CB_Rec_Gen_Energie_Deviation
      << h_IM_CB_interval_Theta_Phi_Energy
      << h_IM_CB_Uncharged_30_Degree_Cut
         //            << h_IM_CB_ZVertex
         //            << h_IM_CB_ZVertex_interval
      << h_IM_CB_ZVertex_interval_30_Degree_Cut
      << h_IM_CB_ZVertex
      << h_Meson_Energy_interval
      << h_Meson_Energy_interval_30_Degree_Cut
      << h_IM_CB_AngleDeviation_Energy
//      << h_IM_CB_AngleDeviation_Photon_Meson_Energy
//      << h_IM_CB_One_high_Photon
      << h_IM_CB_ClusterSize3
      << h_IM_True_Opening_Angle
      << h_IM_Rec_Opening_Angle

         ;

    for( auto h : h_cbs_ClusterSize3) {
        c << h;
    }

    for( auto h : h_cbs_ClusterSize0) {
        c << h;
    }
    c << endc;

}


AUTO_REGISTER_PHYSICS(scratch_sobotzik_Pi0Calib)

