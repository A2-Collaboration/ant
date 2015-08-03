#include "PID_PhiAngle.h"
#include "expconfig/detectors/PID.h"
#include "CalibrationDataManager.h"
#include "detail/Helpers.h"

#include <limits>
#include <cmath>

using namespace ant;
using namespace ant::calibration;
using namespace std;

PID_PhiAngle::ThePhysics::ThePhysics(const string& name, unsigned nChannels) :
    Physics(name)
{
    const BinSettings pid_channels(nChannels);
    const BinSettings phibins(1000, -M_PI, 3*M_PI);

    pid_cb_phi_corr = HistFac.makeTH2D("CB Phi", "IM [MeV]", "#",
                                       phibins, pid_channels, "pid_cb_phi_corr");
}

void PID_PhiAngle::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed().Candidates();

    // search for events with
    // one cluster in CB, one cluster in PID
    // ignore the matched candidates, since this is what
    // we want to calibrate

    const Cluster* cluster_pid = nullptr;
    double phi_cb = numeric_limits<double>::quiet_NaN();

    for(const auto& cand : cands) {

        auto cl_cb_  = cand->FindCluster(detector_t::CB);

        if(cl_cb_ != nullptr) {
            // found more than one CB cluster
            if(isfinite(phi_cb))
                return;
            phi_cb = cand->Phi();
        }

        auto cl_pid_ = cand->FindCluster(detector_t::PID);

        if(cl_pid_ != nullptr) {
            // found more than one pid cluster
            if(cluster_pid != nullptr)
                return;
            cluster_pid = cl_pid_;
        }
    }

    if(!isfinite(phi_cb) || cluster_pid == nullptr)
        return;

    pid_cb_phi_corr->Fill(phi_cb,        cluster_pid->CentralElement);
    pid_cb_phi_corr->Fill(phi_cb+2*M_PI, cluster_pid->CentralElement);

}

void PID_PhiAngle::ThePhysics::Finish()
{

}

void PID_PhiAngle::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << pid_cb_phi_corr << endc;
}

PID_PhiAngle::PID_PhiAngle(shared_ptr<expconfig::detector::PID> pid) :
    Module("PID_PhiAngle"),
    pid_detector(pid)
{
}

PID_PhiAngle::~PID_PhiAngle()
{
}

std::unique_ptr<Physics> PID_PhiAngle::GetPhysicsModule() {
    return std_ext::make_unique<ThePhysics>(GetName(), pid_detector->GetNChannels());
}


std::vector<std::list<TID> > ant::calibration::PID_PhiAngle::GetChangePoints() const
{
    return {calibrationManager->GetChangePoints(GetName())};
}

void ant::calibration::PID_PhiAngle::Update(size_t, const TID& id)
{
    TCalibrationData cdata;
    if(!calibrationManager->GetData(GetName(), id, cdata))
        return;
    if(cdata.Data.size() != 1)
        return;
    const TKeyValue<double>& kv = cdata.Data.front();
    pid_detector->SetPhiOffset(kv.Value);
}
