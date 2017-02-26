#include "ProductionDataBase.h"

using namespace ant;
using namespace ant::mc::data;


static ProductionDataBase::XSections_t::value_type gp_pRho =
{ ParticleTypeTreeDatabase::Channel::gp_pRho,
  ProductionDataBase::MakeInterPolator({
// Eur. Phys. J. A 23, 317–344 (2005)
//
{1080,   0 },
{1100,   6 },
{1150,   12 },
{1250,   16 },
{1400,   20 },
{1600,   25 },
{1800,   22 },
{2000,   20 },
{2200,   18 },
{2400,   16 },
{2600,   16 }
  })};

