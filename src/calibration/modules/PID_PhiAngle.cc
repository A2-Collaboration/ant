#include "PID_PhiAngle.h"
#include "expconfig/detectors/PID.h"
#include "CalibrationDataManager.h"

using namespace ant;
using namespace ant::calibration;
using namespace std;

void PID_PhiAngle::ThePhysics::ProcessEvent(const Event& event)
{

}

void PID_PhiAngle::ThePhysics::Finish()
{

}

void PID_PhiAngle::ThePhysics::ShowResult()
{

}

PID_PhiAngle::PID_PhiAngle(shared_ptr<expconfig::detector::PID> pid) :
    Module("PID_PhiAngle"),
    pid_detector(pid)
{
}

PID_PhiAngle::~PID_PhiAngle()
{
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

}