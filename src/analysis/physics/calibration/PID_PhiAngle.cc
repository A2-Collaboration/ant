#include "PID_PhiAngle.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

PID_PhiAngle::PID_PhiAngle(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    theta_range(40.0*TMath::DegToRad(), 140*TMath::DegToRad())
{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);
    const BinSettings pid_channels(detector->GetNChannels());
    const BinSettings phibins(1000, -180, 3*180);

    pid_cb_phi_corr = HistFac.makeTH2D("CB/PID Cluster/Channel Correlation", "CB Cluster Phi / degree", "#",
                                       phibins, pid_channels, "pid_cb_phi_corr");
}


void PID_PhiAngle::ProcessEvent(const TEvent& event, manager_t&)
{

    // search for events with
    // one cluster in CB, one cluster in PID

    TClusterPtr cluster_pid;
    TClusterPtr cluster_cb;
    for(const TClusterPtr& cl : event.Reconstructed->Clusters) {
        if(cl->DetectorType == Detector_t::Type_t::PID) {
            if(cluster_pid)
                return;
            cluster_pid = cl;
        }
        else if(cl->DetectorType == Detector_t::Type_t::CB) {
            if(cluster_cb)
                return;
            cluster_cb = cl;
        }
    }

    if(cluster_pid && cluster_cb) {
        const double phi_cb_degrees = std_ext::radian_to_degree(cluster_cb->Position.Phi());

        pid_cb_phi_corr->Fill(phi_cb_degrees,     cluster_pid->CentralElement);
        pid_cb_phi_corr->Fill(phi_cb_degrees+360, cluster_pid->CentralElement);
    }
}

void PID_PhiAngle::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << pid_cb_phi_corr << endc;
}

AUTO_REGISTER_PHYSICS(PID_PhiAngle)

