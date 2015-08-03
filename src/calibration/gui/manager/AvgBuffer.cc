#include "AvgBuffer.h"

#include "TH1D.h"
#include "TH2D.h"
#include "base/interval.h"
#include "tree/TDataRecord.h"

using namespace ant::calibration::gui;

template class AvgBuffer<TH1D,ant::interval<ant::TID>>;
template class AvgBuffer<TH2D,ant::interval<ant::TID>>;
