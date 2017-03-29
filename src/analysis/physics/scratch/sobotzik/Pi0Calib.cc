#include "Pi0Calib.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/Logger.h"
#include "TH1D.h"
#include "TTree.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

scratch_sobotzik_Pi0Calib::scratch_sobotzik_Pi0Calib(const string& name, OptionsPtr opts)
    : Physics(name, opts)
{
    const auto& caloEnergyWindow = opts->Get<interval<double>>("CaloEnergyWindow", {-std_ext::inf, std_ext::inf});
    hists.emplace_back(hist_t{HistFac, {2,2}, hist_t::any,caloEnergyWindow});
//    hists.emplace_back(hist_t{HistFac, {3,3}, hist_t::any},caloEnergyWindow});
//    hists.emplace_back(hist_t{HistFac, {5,5}, hist_t::any},caloEnergyWindow});
//    hists.emplace_back(hist_t{HistFac, {6,6}, hist_t::any},caloEnergyWindow});

    const BinSettings bins_Mult(10);

    h_Mult_All    = HistFac.makeTH1D("Multiplicity: All",   "n Clusters/Event","", bins_Mult,"n_All");
    h_Mult_CB     = HistFac.makeTH1D("Multiplicity: CB",    "n Clusters/Event","", bins_Mult,"n_CB");

    h_Mult_TAPS   = HistFac.makeTH1D("Multiplicity: TAPS",  "n Clusters/Event","", bins_Mult,"n_TAPS");
}

const scratch_sobotzik_Pi0Calib::hist_t::range_t scratch_sobotzik_Pi0Calib::hist_t::any = {0, numeric_limits<int>::max()};

scratch_sobotzik_Pi0Calib::hist_t::hist_t(const HistogramFactory& HistFac,
                                 const range_t& cb, const range_t& taps
                                 ,const interval<double>& caloEnergy_window) :
    n_CB(cb), n_TAPS(taps),CaloEnergy_Window(caloEnergy_window)
{
    auto to_string = [] (const range_t& r) {
        if(r == any)
            return string("any");
        if(r.Start() == r.Stop())
            return std::to_string(r.Start());
        return std::to_string(r.Start())+std::to_string(r.Stop());
    };

    prefix = std_ext::formatter()
                   << "h_"
                   << to_string(n_CB)   << "CB_"
                   << to_string(n_TAPS) << "TAPS";

    HistogramFactory histFac(prefix, HistFac, prefix);
    const BinSettings bins_IM   (400, 0, 1100); // MeV
    const BinSettings bins_angle(180, 0,    180); // degrees
    const BinSettings bins_timing(200,-100,100); // ns
    const BinSettings bins_energy(200, 0, 1000);

    h_IM_All   = histFac.makeTH1D("IM: All",  "IM / MeV","",bins_IM,"IM_All");

    h_IM_CB_all             = histFac.makeTH2D("IM: CB",   "IM / MeV","E [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_All");
    h_IM_CB_Uncharged_No_Cut             = histFac.makeTH2D("IM: CB",   "IM / MeV","E [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Uncharged");
    h_IM_CB_interval        = histFac.makeTH2D("IM: CB",   "IM / MeV","E [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Interval");
    h_IM_CB_interval_Uncharged_No_Cut        = histFac.makeTH2D("IM: CB",   "IM / MeV","E [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Interval_No_Cut");
    h_IM_CB_interval_Uncharged_30_Degree_Cut        = histFac.makeTH2D("IM: CB",   "IM / MeV","E [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Interval_30_Degree_Cut");

    h_IM_CB_Uncharged_30_Degree_Cut    = histFac.makeTH2D("IM: CB",   "IM / MeV","E [MeV]",bins_IM,BinSettings(32,0,800),"IM_CB_Uncharged_30_Degree_Cut");

    h_IM_CB_Angle_Energy    = histFac.makeTH2D("IM: Angle",   "Angle / Degrees","E [MeV]",bins_angle,BinSettings(32,0,800),"IM_CB_Angle");
    h_IM_CB_interval_Theta_Phi_Energy= histFac.makeTH3D("IM:CB","Polar angle Theta / Degree","Azimut angle Phi / Degree","Energz of the Photons E[MeV]", bins_angle,BinSettings(360,-180,180) ,BinSettings(32,0,800),"IM_CB_Interval_Theta_Phi");

    h_IM_CB_ZVertex         = histFac.makeTH3D("IM: CB",   "IM / MeV","E [MeV]","Z-Vertex [cm]",bins_IM,BinSettings(32,0,800),BinSettings(10,-5,5),"IM_CB_ZVertex");
    h_IM_CB_ZVertex_interval         = histFac.makeTH3D("IM: CB",   "IM / MeV","E [MeV]","Z-Vertex [cm]",bins_IM,BinSettings(32,0,800),BinSettings(10,-5,5),"IM_CB_ZVertex_interval");
    h_IM_CB_ZVertex_interval_30_Degree_Cut         = histFac.makeTH3D("IM: CB",   "IM / MeV","E [MeV]","Z-Vertex [cm]",bins_IM,BinSettings(32,0,800),BinSettings(10,-5,5),"IM_CB_ZVertex_interval_30_Degree_Cut");


    h_IM_CB_corr    = histFac.makeTH1D("IM: CB corr",   "IM / MeV","",bins_IM,"IM_CB_corr");
    h_IM_TAPS  = histFac.makeTH1D("IM: TAPS", "IM / MeV","",bins_IM,"IM_TAPS");

    h_Angle_CB   = histFac.makeTH1D("Angle: CB",   "angle [#circ]","",bins_angle,"Angle_CB");
    h_Angle_TAPS = histFac.makeTH1D("Angle: TAPS", "angle [#circ]","",bins_angle,"Angle_TAPS");



    h_ClusterHitTiming_CB   = histFac.makeTH2D("ClusterHitTiming: CB",   "Energy","t / ns",bins_energy,bins_timing,"ClusterHitTiming_CB");
    h_ClusterHitTiming_TAPS = histFac.makeTH2D("ClusterHitTiming: TAPS", "Energy","t / ns",bins_energy,bins_timing,"ClusterHitTiming_TAPS");

}

void scratch_sobotzik_Pi0Calib::hist_t::Fill(const TCandidatePtrList& c_CB, const TCandidatePtrList& c_TAPS, const double zVertex) const
{
    if(!n_CB.Contains(c_CB.size()))
        return;
    if(!n_TAPS.Contains(c_TAPS.size()))
        return;

    //LOG(INFO) << c_CB.at(0)->Theta;
//    double angleedge = 30;
//    if  (c_CB.at(0)->Theta <(angleedge * 2 * 3.141 /360) ||c_CB.at(0)->Theta >180 - (angleedge * 2 * 3.141 /360))
//    {
//        return;
//    }
//    else
//    {
//        if (c_CB.at(1)->Theta <(angleedge * 2 * 3.141 /360)|| c_CB.at(1)->Theta > 180 - (angleedge * 2 * 3.141 /360))
//        {
//            return;
//        }
//    }

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
    double angleedge = 30;
    const auto& sum_CB = sum_as_photons(c_CB);
    const auto& sum_TAPS = sum_as_photons(c_TAPS);
    const auto& angle_CB = min_angle(c_CB);

    h_IM_All->Fill((sum_CB+sum_TAPS).M());

    const auto bin1 = h_IM_CB_interval->GetYaxis()->FindBin(c_CB.at(0)->CaloEnergy);
    const auto bin2 = h_IM_CB_interval->GetYaxis()->FindBin(c_CB.at(1)->CaloEnergy);
    if(bin1==bin2) {
        if(sum_CB.M()>1.0) {
            h_IM_CB_interval->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
            h_IM_CB_ZVertex_interval->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,zVertex);
            if((c_CB.at(0)->VetoEnergy == 0 )&&(c_CB.at(1)->VetoEnergy == 0))
            {
                h_IM_CB_interval_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
                h_IM_CB_interval_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy);


                h_IM_CB_interval_Theta_Phi_Energy->Fill(c_CB.at(0)->Theta / (2 * 3.141) *360,c_CB.at(0)->Phi / (2 * 3.141) *360,c_CB.at(0)->CaloEnergy);
                h_IM_CB_interval_Theta_Phi_Energy->Fill(c_CB.at(1)->Theta / (2 * 3.141) *360,c_CB.at(1)->Phi / (2 * 3.141) *360 ,c_CB.at(1)->CaloEnergy);

                if(     (c_CB.at(0)->Theta >(angleedge * 2 * 3.141 /360) &&
                         c_CB.at(0)->Theta <180 - (angleedge * 2 * 3.141 /360))
                        &&
                        (c_CB.at(1)->Theta >(angleedge * 2 * 3.141 /360) &&
                         c_CB.at(1)->Theta <180 - (angleedge * 2 * 3.141 /360)))
                {
                    h_IM_CB_interval_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(0)->CaloEnergy);
                    h_IM_CB_interval_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(1)->CaloEnergy);
                    h_IM_CB_ZVertex_interval_30_Degree_Cut->Fill (sum_CB.M(),c_CB.at(0)->CaloEnergy,zVertex);
                }
            }

        }
    }

//    zVertex..
//    LOG(INFO)<< c_CB.at(0)->Phi;


    if(sum_CB.M()>1.0)
    {
        h_IM_CB_all->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
        h_IM_CB_all->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy);
        h_IM_CB_Angle_Energy->Fill( angle_CB,c_CB.at(0)->CaloEnergy);
        h_IM_CB_Angle_Energy->Fill( angle_CB,c_CB.at(1)->CaloEnergy);
        h_IM_CB_ZVertex->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy,zVertex);

        if((c_CB.at(0)->VetoEnergy == 0 )&&(c_CB.at(1)->VetoEnergy == 0))
        {
            h_IM_CB_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(0)->CaloEnergy);
            h_IM_CB_Uncharged_No_Cut->Fill(sum_CB.M(),c_CB.at(1)->CaloEnergy);


            if(     (c_CB.at(0)->Theta >(angleedge * 2 * 3.141 /360) &&
                     c_CB.at(0)->Theta <180 - (angleedge * 2 * 3.141 /360))
                    &&
                    (c_CB.at(1)->Theta >(angleedge * 2 * 3.141 /360) &&
                     c_CB.at(1)->Theta <180 - (angleedge * 2 * 3.141 /360)))
            {
                h_IM_CB_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(0)->CaloEnergy);
                h_IM_CB_Uncharged_30_Degree_Cut->Fill( sum_CB.M(),c_CB.at(1)->CaloEnergy);
            }
        }


    }



    h_IM_CB_corr->Fill(sum_as_corr_photons(c_CB).M());
    h_IM_TAPS->Fill(sum_TAPS.M());

    h_Angle_CB->Fill(min_angle(c_CB));
    h_Angle_TAPS->Fill(min_angle(c_TAPS));

    fill_timing(c_CB, h_ClusterHitTiming_CB);
    fill_timing(c_TAPS, h_ClusterHitTiming_TAPS);
}

void scratch_sobotzik_Pi0Calib::hist_t::ShowResult() const
{
    canvas(prefix)
//            << h_Angle_CB
//            << h_Angle_TAPS
//            << h_IM_All
            << drawoption("colz")
            << h_IM_CB_all
            << h_IM_CB_interval
            << h_IM_CB_interval_Uncharged_No_Cut
            << h_IM_CB_interval_Uncharged_30_Degree_Cut
            << h_IM_CB_Uncharged_No_Cut
            << h_IM_CB_Angle_Energy
            << h_IM_CB_interval_Theta_Phi_Energy
            << h_IM_CB_Uncharged_30_Degree_Cut
            << h_IM_CB_ZVertex
            << h_IM_CB_ZVertex_interval
            << h_IM_CB_ZVertex_interval_30_Degree_Cut
            << endc;
}


scratch_sobotzik_Pi0Calib::~scratch_sobotzik_Pi0Calib()
{

}

void scratch_sobotzik_Pi0Calib::ProcessEvent(const TEvent& event, manager_t&)
{

    auto ptree = event.MCTrue().ParticleTree;
    if(ptree) {
        auto typetree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g);
        if(!ptree->IsEqual(typetree, utils::ParticleTools::MatchByParticleName))
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

    for(auto& h : hists)
        h.Fill(c_CB, c_TAPS, event.MCTrue().Target.Vertex.z);

    h_Mult_All->Fill(event.Reconstructed().Candidates.size());
    h_Mult_CB->Fill(c_CB.size());
    h_Mult_TAPS->Fill(c_TAPS.size());
}

void scratch_sobotzik_Pi0Calib::ShowResult()
{
    for(const auto& h : hists) {
        h.ShowResult();
    }

    canvas(GetName())
            << h_Mult_All
            << h_Mult_CB
            << h_Mult_TAPS
            << endc;
}


AUTO_REGISTER_PHYSICS(scratch_sobotzik_Pi0Calib)

