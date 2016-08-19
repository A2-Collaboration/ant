#include "ExtractShowerDepth.h"

#include "plot/root_draw.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

ExtractShowerDepth::ExtractShowerDepth(const string& name, OptionsPtr opts) :
    Physics(name, opts)
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

    t.TrueZVertex = mctrue.Target.Vertex.z;
    t.TrueTheta = std_ext::radian_to_degree(true_particle->Theta());
    t.TrueEk = true_particle->Ek();
    t.Theta = std_ext::radian_to_degree(calocluster->Position.Theta());
    t.Ek = calocluster->Energy;
    t.Tree->Fill();

}

void ExtractShowerDepth::ShowResult()
{
    canvas(GetName()) << steps
                      << drawoption("colz")
                      << TTree_drawable(t.Tree, "Ek:TrueEk >> h1__(1000,0,1600,1000,0,1600)","")
                      << TTree_drawable(t.Tree, "Theta:TrueTheta >> h2__(180,0,180,180,0,180)","")
                      << endc;
    canvas(GetName()+" Z Vertex")
            << drawoption("colz")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h1(180,0,180,40,-10,10)","")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h2(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h3(180,0,180,40,-10,10)","-2<TrueZVertex && TrueZVertex<-1")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h4(180,0,180,40,-10,10)","-1<TrueZVertex && TrueZVertex<+0")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h5(180,0,180,40,-10,10)","+0<TrueZVertex && TrueZVertex<+1")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h6(180,0,180,40,-10,10)","+1<TrueZVertex && TrueZVertex<+2")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h7(180,0,180,40,-10,10)","+2<TrueZVertex && TrueZVertex<+3")
            << endc;
    canvas(GetName()+" ZVertex Energy")
            << drawoption("colz")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h1_(180,0,180,40,-10,10)","")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h2_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 000<TrueEk && TrueEk<100")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h3_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 100<TrueEk && TrueEk<200")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h4_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 200<TrueEk && TrueEk<300")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h5_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 300<TrueEk && TrueEk<400")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h6_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 400<TrueEk && TrueEk<500")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h7_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 500<TrueEk && TrueEk<600")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h8_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 600<TrueEk && TrueEk<700")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> h9_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 700<TrueEk && TrueEk<800")
            << TTree_drawable(t.Tree, "(Theta-TrueTheta):TrueTheta >> ha_(180,0,180,40,-10,10)","-3<TrueZVertex && TrueZVertex<-2 && 800<TrueEk && TrueEk<900")
            << endc;
}

AUTO_REGISTER_PHYSICS(ExtractShowerDepth)
