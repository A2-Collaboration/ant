#include "PID_PhiAngle.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

PID_PhiAngle::PID_PhiAngle(const string& name, analysis::PhysOptPtr opts) :
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
    const auto& cands = event.Reconstructed->Candidates;

    // search for events with
    // one cluster in CB, one cluster in PID
    // ignore the matched candidates, since this is what
    // we want to calibrate

    TClusterPtr cluster_pid = nullptr;
    double phi_cb = numeric_limits<double>::quiet_NaN();

    for(const auto& cand : cands) {

        if(!theta_range.Contains(cand->Theta))
            continue;

        auto cl_cb_  = cand->FindFirstCluster(Detector_t::Type_t::CB);

        if(cl_cb_ != nullptr) {
            // found more than one CB cluster
            if(isfinite(phi_cb))
                return;
            phi_cb = cand->Phi;
        }

        auto cl_pid_ = cand->FindFirstCluster(Detector_t::Type_t::PID);

        if(cl_pid_ != nullptr) {
            // found more than one pid cluster
            if(cluster_pid != nullptr)
                return;
            cluster_pid = cl_pid_;
        }
    }

    /// \todo search all clusters, leave candidates alone
    for(const TClusterPtr& cl : event.Reconstructed->Clusters) {
        if(cl->DetectorType != Detector_t::Type_t::PID)
            continue;
        if(!isfinite(cl->Energy) || !isfinite(cl->Time))
            continue;
        // found more than one PID cluster
        if(cluster_pid != nullptr)
            return;
        cluster_pid = cl;
    }

    if(!isfinite(phi_cb) || cluster_pid == nullptr)
        return;

    const double phi_cb_degrees = std_ext::radian_to_degree(phi_cb);

    pid_cb_phi_corr->Fill(phi_cb_degrees,     cluster_pid->CentralElement);
    pid_cb_phi_corr->Fill(phi_cb_degrees+360, cluster_pid->CentralElement);
}

void PID_PhiAngle::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << pid_cb_phi_corr << endc;
}

AUTO_REGISTER_PHYSICS(PID_PhiAngle)

