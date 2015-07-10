ant
===

Just another **AN**alysis **T**oolkit `ant`,
based on GoAT tree output from AcquRoot.

It also includes it's own unpacker for acqu data files,
but this is WIP.

## Dependencies
  * C++11
  * cmake
  * doxygen
  * [CERN ROOT5](https://root.cern.ch/)
  * [PLUTO](https://www-hades.gsi.de/?q=pluto)
  * [APLCON++](https://github.com/neiser/APLCON)

## External Components
  * [Easylogging++](http://easylogging.muflihun.com/)
  * [Catch] (https://github.com/philsquared/Catch) (framework for unit-tests, TDD)

## GoAT File Reader

Almost implemented.

#### Tree Manager
TODOs:
 * remove tree is module that requested it fails to init. Right now only the module is removed, leaving the tree active
 * keep track of associated branches -> Allow modules to share trees
 * unset branch addresses on module destruction? -> Possible (probable) segfault here, if module terminates after setting a branch and tree wants to read into it
 * unset unused branches? -> speed?

#### Particles
Currently negelecting energy and angle inforamtion from particle trees
and just using the referenced track. This info is a dublicate of the
track data anyway. Ohterwise no idea what to do with it.

#### Detector Hits
Ant data structure missing.

#### Pluto
Works, but no decay tree inforamtion. PParticle::GetDauther() always returns nullptr.


## Brainstorm stuff
ant
 base
   - Logger
   - printable
   - stl extensions

  analysis
   - datastructures
   - plot
   - utils
   - input
     - goat
     - unpacker
   - physics classes

   unpacker
    - raw file reader
    - unpacker cfg
    - logical channel mapping

   reconstruct
    - apply "calibration"
    - clustering
    - track matching
   
   Exp cfg
    - ADC <-> log. channel cfg
    - element positions, properties...
    - broken channels

   Calibrations


