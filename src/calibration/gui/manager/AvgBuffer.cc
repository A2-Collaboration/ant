#include "AvgBuffer.h"

#include "TH1.h"
#include "base/interval.h"
#include "tree/TDataRecord.h"

using namespace ant::calibration::gui;

template class AvgBuffer<TH1,ant::interval<ant::TID>>;
