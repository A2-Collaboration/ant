#include "Reconstruct.h"

#include "HitMatching.h"
#include "Clustering.h"
#include "ApplyCalibrations.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;

Reconstruct::Reconstruct(const THeaderInfo &headerInfo)
{
  auto config = ExpConfig::Reconstruct::Get(headerInfo);
}

unique_ptr<TEvent> Reconstruct::DoReconstruct(unique_ptr<TDetectorRead>& read)
{

}
