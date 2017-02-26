#include "ProductionDataBase.h"

using namespace ant;
using namespace ant::mc::data;


static ProductionDataBase::XSections_t::value_type  gp_pOmega =
{ ParticleTypeTreeDatabase::Channel::gp_pOmega,
  ProductionDataBase::MakeInterPolator({
    { 1110, 0.000 },
    { 1112, 0.999 },
    { 1138, 4.017 },
    { 1163, 5.989 },
    { 1187, 7.213 },
    { 1225, 7.948 },
    { 1275, 8.508 },
    { 1325, 7.755 },
    { 1375, 7.185 },
    { 1424, 6.724 },
    { 1474, 6.307 },
    { 1525, 6.498 },
    { 1575, 6.178 },
    { 1625, 6.309 },
    { 1675, 6.752 },
    { 1750, 6.361 },
    { 1850, 6.677 },
    { 1950, 6.790 },
    { 2050, 6.498 },
    { 2150, 6.145 },
    { 2250, 5.636 },
    { 2349, 5.509 },
    { 2450, 4.850 },
    { 2550, 4.784 }
  })};

