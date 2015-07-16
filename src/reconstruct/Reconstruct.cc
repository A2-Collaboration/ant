#include "Reconstruct.h"

#include "HitMatching.h"
#include "Clustering.h"
#include "ApplyCalibrations.h"

#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"

using namespace std;
using namespace ant;

Reconstruct::Reconstruct(const THeaderInfo &headerInfo)
{
  auto config = ExpConfig::Reconstruct::Get(headerInfo);
  calibrations = config->GetCalibrations();
  updateables = config->GetUpdateables();
}

unique_ptr<TEvent> Reconstruct::DoReconstruct(TDetectorRead& read)
{

}
