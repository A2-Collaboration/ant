[![Build status master branch](https://travis-ci.org/A2-Collaboration-dev/ant.svg?branch=master)](https://travis-ci.org/A2-Collaboration-dev/ant)
[![Codecov status master branch](https://codecov.io/github/A2-Collaboration-dev/ant/coverage.svg?branch=master)](https://codecov.io/github/A2-Collaboration-dev/ant?branch=master)


ant
===

Just another **AN**alysis **T**oolkit `ant`, which can read from many
input sources with minimal user intervention and let's you create
Physics analysis within minutes.

Please see also the automatically generated
[Doxygen pages](http://a2-collaboration-dev.github.io/ant/).



## Dependencies
  * C++11 (gcc >4.8.2 should work)
  * cmake
  * doxygen (optional)
  * [CERN ROOT5](https://root.cern.ch/)
  * [PLUTO](https://www-hades.gsi.de/?q=pluto)
  * [APLCON++](https://github.com/neiser/APLCON)

# Contributing

Please read the following sections if you want to contribute to this
project. 

## Coding Style
  * Indentation: 4 spaces, no tabs anywhere
  * Ant-Codesytle defined in [doc/Ant-Codestyle.xml](doc/Ant-Codestyle.xml) (import to QtCreator)
  * `#include` statements:
    * header file for this .cc file first
    * then grouped by ant, ROOT, STL, others
    * one line between groups
    * each group ordered alphabetically

## Development Status
  * [x] Unpacker for Acqu Mk2 data
  * [x] Unpacker for a2geant data
  * [ ] [Calibration modules](src/calibration/modules):
    * [x] CB Energy/Timewalk
    * [x] TAPS Energy
    * [x] Time offsets for PID/CB/TAPS/Tagger
    * [ ] Complete calibration cycle on blaster
  * [ ] [Experiment configuration](src/expconfig/setups)
    * [x] EPT 2014 beamtimes
    * [ ] Any normal tagger beamtime
    * [ ] Wire chamber
  * [x] Reconstruct
    * [x] Apply calibration factors
    * [x] Update calibration factors
    * [x] Clustering
    * [x] Candidate builder (Veto/Calorimeter matching)
  * [ ] Physics
    * [x] Data structure for events
    * [x] Input readers
    * [ ] Slowcontrol handling
    * [ ] Kinematic fitting on measured data


## Detector Type Mapping

| Ant Reconstruct  | Ant Analysis  | Goat |
|------------------|---------------|------|
| Trigger          | -             | -    |
| Tagger           | -             | -    |
| TaggerMicroscope | -             | -    |
| EPT              | -             | -    |
| Moeller          | -             | -    |
| PairSpec         | -             | -    |
| CB               | CB            | NaI  |
| PID              | PID           | PID  |
| MWPC0            | MWPC0         | MWPC |
| MWPC1            | MWPC1         |      |
| TAPS             | TAPS          | BaF2 |
|                  |               | PbWO4|
| TAPSVeto         | TAPSVeto      | Veto |
| Cherenkov        | Cherenkov     | -    |

## External Components

Have a look at those very nice projects, which are used here:

  * [Easylogging++](http://easylogging.muflihun.com/)
  * [Catch](https://github.com/philsquared/Catch) framework for unit-tests, test-driven development. See [the test/ subdirectory](test/).
  * [TCLAP - Templatized C++ Command Line Parser](http://tclap.sourceforge.net)
