# ant::analysis

## Unit System
  * Energy: MeV
  * Angle: radians
  * Time: ns
  * Distance: cm? (not used atm. make conform with the rest of the project)

## Data Structure Philosophy
The central data structure in the analysis is the ant::Event. It contains all inforamtion required to analyse on an event basis:
Tracks, Particles, TaggerHits, Trigger Infos, each for recontructed and mctrue.
Eeach of the above groups is represented as a list/vector of shared pointers to the data classes. This solves the owenship problem in an elegant way.
The user does not have to take care about deleting and object lifetime.

## TODO
  * Move all analysis stuff to ant::analysis namespace.
  * ant input (to read stuff from the unpacker/reconstruct chain)
    * find useful conversion
    * maybe redesing stuff to be closer to the rest of ant / more common
      * use same Detector_t ?
      * keep using shared_ptr for everything or move to vector<vector<>> based style as the rest? (and maybe avoid copying every data object to a new format)
  * Investigate memory leak due to circular ref in Particle relationships




