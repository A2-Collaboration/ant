#include "ProductionDataBase.h"

using namespace ant;
using namespace ant::mc::data;


// from PHYSICAL REVIEW C 88, 044601 (2013)
// Measurement of the γ p → K 0 Sigma + reaction with the Crystal Ball/TAPS detectors at the Mainz Microtron
// for higher energies: SAPHIR @ CB/ELSA (https://inspirehep.net/record/680746)
// Eur.Phys.J. A24 (2005) 275-286
// Measurement of the reaction gamma p ---> K0 Sigma+ at photon energies up to 2.6-GeV

static ProductionDataBase::XSections_t::value_type gp_SigmaPlusK0S =
{ ParticleTypeTreeDatabase::Channel::gp_SigmaPlusK0S,
  ProductionDataBase::MakeInterPolator({

    { 1075.0,    0.0},    // note: changed manually!!! was -0.1 in original data!
    { 1100.0,   .135},
    { 1125.0,   .192},
    { 1150.0,   .217},
    { 1175.0,   .233},
    { 1200.0,   .280},
    { 1225.0,   .325},
    { 1250.0,   .385},
    { 1275.0,   .448},
    { 1300.0,   .512},
    { 1325.0,   .580},
    { 1350.0,   .633},
    { 1375.0,   .675},
    { 1400.0,   .723},
    { 1425.0,   .695},
    { 1500.0,   .588},  // values from SAPHIR paper, bin center chosen as energy value
    { 1600.0,   .575},
    { 1700.0,   .466},
    { 1800.0,   .391},
    { 1950.0,   .238},
    { 2150.0,   .201},
    { 2350.0,   .125},
    { 2550.0,   .086}


  })};

