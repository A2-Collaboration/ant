#include "ExtractShowerDepth.h"

#include "plot/root_draw.h"
#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/TAPS.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

ExtractShowerDepth::ExtractShowerDepth(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    MaxTheta(opts->Get<double>("MaxTheta", 180.0)),
    CB_param1(opts->Get<double>("CB_param1", 3.0))
{

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(10),"steps");
    t.CreateBranches(HistFac.makeTTree("tree"));
}

void ExtractShowerDepth::ProcessEvent(const TEvent& event, manager_t&)
{
    const TEventData& recon = event.Reconstructed();
    const TEventData& mctrue = event.MCTrue();

    steps->Fill("Seen",1);

    if(mctrue.Particles.GetAll().size() != 1)
        return;
    steps->Fill("MC True nParticles==1",1);

    auto& true_particle = mctrue.Particles.GetAll().front();

    TClusterPtr calocluster;

    for(auto& cluster : recon.Clusters.get_iter()) {
        if(cluster->DetectorType & Detector_t::Any_t::Calo)
            calocluster = cluster;
    }
    if(!calocluster)
        return;
    steps->Fill("nCaloClusters==1",1);

    auto theta = calocluster->Position.Theta();
    auto theta_true = true_particle->Theta();
    auto Ek = calocluster->Energy;
    auto z_vertex = mctrue.Target.Vertex.z;

    t.TrueZVertex = z_vertex;
    t.TrueTheta = std_ext::radian_to_degree(theta_true);
    t.TrueEk = true_particle->Ek();
    t.Theta = std_ext::radian_to_degree(theta);
    t.Ek = Ek;

    // copied from Fitter::FitParticle::GetVector
    t.ThetaCorr = std_ext::NaN;
    t.ShowerDepth = std_ext::NaN;
    t.RadiationLength = std_ext::NaN;
    if(calocluster->DetectorType == Detector_t::Type_t::CB) {
        static auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
        auto elem = cb->GetClusterElement(calocluster->CentralElement);

        t.ShowerDepth = z_vertex / (std::cos(theta) - std::cos(theta_true)) - cb->GetInnerRadius();
        t.RadiationLength = t.ShowerDepth / (elem->RadiationLength*std::log2(Ek/elem->CriticalE));

        const auto R  = cb->GetInnerRadius() + 1.0*elem->RadiationLength*std::log2(Ek/elem->CriticalE)/std::pow(std::sin(theta),CB_param1);
        t.ThetaCorr = std::acos(( R*std::cos(theta) - z_vertex) / R );
    }
    else if(calocluster->DetectorType == Detector_t::Type_t::TAPS) {
        static auto taps = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
        auto elem = taps->GetClusterElement(calocluster->CentralElement);
        const auto Z  =  taps->GetZPosition() + elem->RadiationLength*std::log2(Ek/elem->CriticalE);
        t.ThetaCorr = std::atan( Z*std::tan(theta) / (Z - z_vertex));
    }
    t.ThetaCorr = std_ext::radian_to_degree(t.ThetaCorr());

    t.Tree->Fill();

}

void ExtractShowerDepth::ShowResult()
{
    const string bins_theta = "180,0,"+to_string(MaxTheta);
    canvas("Overview") << steps
                       << drawoption("colz")
                       << TTree_drawable(t.Tree, "Ek:TrueEk >> (1000,0,1600,1000,0,1600)","")
                       << TTree_drawable(t.Tree, "(Ek-TrueEk)/TrueEk >> (200,-0.2,0.2)","")
                       << TTree_drawable(t.Tree, "(1.0/Ek-1.0/TrueEk)/(1.0/TrueEk) >> (200,-0.2,0.2)","")
                       << TTree_drawable(t.Tree, "Theta:TrueTheta >>("+bins_theta+","+bins_theta+")","")
                       << TTree_drawable(t.Tree, "ThetaCorr:TrueTheta >>("+bins_theta+","+bins_theta+")","")
                       << endc;
    const string bins_theta_op = ">>("+bins_theta+",";
    {
        string formula = "(Theta-TrueTheta):TrueTheta";
        string bins = bins_theta_op+"40,-10,10)";
        canvas("Theta")
                << drawoption("colz")
                << TTree_drawable(t.Tree, formula+bins,"")
                << TTree_drawable(t.Tree, formula+bins,"-5<TrueZVertex && TrueZVertex<-4")
                << TTree_drawable(t.Tree, formula+bins,"-4<TrueZVertex && TrueZVertex<-3")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2")
                << TTree_drawable(t.Tree, formula+bins,"-2<TrueZVertex && TrueZVertex<-1")
                << TTree_drawable(t.Tree, formula+bins,"-1<TrueZVertex && TrueZVertex<+0")
                << TTree_drawable(t.Tree, formula+bins,"+0<TrueZVertex && TrueZVertex<+1")
                << TTree_drawable(t.Tree, formula+bins,"+1<TrueZVertex && TrueZVertex<+2")
                << TTree_drawable(t.Tree, formula+bins,"+2<TrueZVertex && TrueZVertex<+3")
                << TTree_drawable(t.Tree, formula+bins,"+3<TrueZVertex && TrueZVertex<+4")
                << TTree_drawable(t.Tree, formula+bins,"+4<TrueZVertex && TrueZVertex<+5")
                << endc;
    }
    {
        string formula = "(ThetaCorr-TrueTheta):TrueTheta";
        string bins = bins_theta_op+"40,-10,10)";
        canvas("ThetaCorr "+to_string(CB_param1))
                << drawoption("colz")
                << TTree_drawable(t.Tree, formula+bins,"")
                << TTree_drawable(t.Tree, formula+bins,"-5<TrueZVertex && TrueZVertex<-4")
                << TTree_drawable(t.Tree, formula+bins,"-4<TrueZVertex && TrueZVertex<-3")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2")
                << TTree_drawable(t.Tree, formula+bins,"-2<TrueZVertex && TrueZVertex<-1")
                << TTree_drawable(t.Tree, formula+bins,"-1<TrueZVertex && TrueZVertex<+0")
                << TTree_drawable(t.Tree, formula+bins,"+0<TrueZVertex && TrueZVertex<+1")
                << TTree_drawable(t.Tree, formula+bins,"+1<TrueZVertex && TrueZVertex<+2")
                << TTree_drawable(t.Tree, formula+bins,"+2<TrueZVertex && TrueZVertex<+3")
                << TTree_drawable(t.Tree, formula+bins,"+3<TrueZVertex && TrueZVertex<+4")
                << TTree_drawable(t.Tree, formula+bins,"+4<TrueZVertex && TrueZVertex<+5")
                << endc;
    }
    {
        string formula = "ShowerDepth:TrueTheta";
        string bins = bins_theta_op+"40,0,30)";
        canvas("ShowerDepth")
                << drawoption("colz")
                << TTree_drawable(t.Tree, formula+bins,"")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2")
                << TTree_drawable(t.Tree, formula+bins,"-2<TrueZVertex && TrueZVertex<-1")
                << TTree_drawable(t.Tree, formula+bins,"-1<TrueZVertex && TrueZVertex<+0")
                << TTree_drawable(t.Tree, formula+bins,"+0<TrueZVertex && TrueZVertex<+1")
                << TTree_drawable(t.Tree, formula+bins,"+1<TrueZVertex && TrueZVertex<+2")
                << TTree_drawable(t.Tree, formula+bins,"+2<TrueZVertex && TrueZVertex<+3")
                << endc;
    }
    {
        string formula = "RadiationLength:Theta";
        string bins = bins_theta_op+"40,0,5)";
        canvas("RadiationLength")
                << drawoption("colz")
                << TTree_drawable(t.Tree, formula+bins,"")
                << TTree_drawable(t.Tree, formula+bins,"-5<TrueZVertex && TrueZVertex<-4")
                << TTree_drawable(t.Tree, formula+bins,"-4<TrueZVertex && TrueZVertex<-3")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2")
                << TTree_drawable(t.Tree, formula+bins,"-2<TrueZVertex && TrueZVertex<-1")
                << TTree_drawable(t.Tree, formula+bins,"-1<TrueZVertex && TrueZVertex<+0")
                << TTree_drawable(t.Tree, formula+bins,"+0<TrueZVertex && TrueZVertex<+1")
                << TTree_drawable(t.Tree, formula+bins,"+1<TrueZVertex && TrueZVertex<+2")
                << TTree_drawable(t.Tree, formula+bins,"+2<TrueZVertex && TrueZVertex<+3")
                << TTree_drawable(t.Tree, formula+bins,"+3<TrueZVertex && TrueZVertex<+4")
                << TTree_drawable(t.Tree, formula+bins,"+4<TrueZVertex && TrueZVertex<+5")
                << endc;
    }
    {
        string formula = "(ThetaCorr-TrueTheta):TrueTheta";
        string bins = bins_theta_op+"40,-10,10)";
        canvas("ThetaCorr Energy")
                << drawoption("colz")
                << TTree_drawable(t.Tree, formula+bins,"")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  000<TrueEk && TrueEk<100")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  100<TrueEk && TrueEk<200")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  200<TrueEk && TrueEk<300")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  300<TrueEk && TrueEk<400")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  400<TrueEk && TrueEk<500")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  500<TrueEk && TrueEk<600")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  600<TrueEk && TrueEk<700")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  700<TrueEk && TrueEk<800")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  800<TrueEk && TrueEk<900")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 &&  900<TrueEk && TrueEk<1000")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 && 1000<TrueEk && TrueEk<1100")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 && 1100<TrueEk && TrueEk<1200")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 && 1200<TrueEk && TrueEk<1300")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 && 1300<TrueEk && TrueEk<1400")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 && 1400<TrueEk && TrueEk<1500")
                << TTree_drawable(t.Tree, formula+bins,"-3<TrueZVertex && TrueZVertex<-2 && 1500<TrueEk && TrueEk<1600")
                << endc;
    }
}

AUTO_REGISTER_PHYSICS(ExtractShowerDepth)
