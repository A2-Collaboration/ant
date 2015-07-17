#ifdef __CINT__
// some ROOTCINT defaults
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// the relevant classes
#pragma link C++ class ant::TID+;
#pragma link C++ class ant::TDataRecord+;
#pragma link C++ class ant::TDetectorReadHit+;
#pragma link C++ class ant::TDetectorRead+;
#pragma link C++ class ant::TEvent+;
#pragma link C++ class ant::TClusterHitDatum+;
#pragma link C++ class ant::TClusterHit+;
#pragma link C++ class ant::TCluster+;
#pragma link C++ class ant::TTrack+;
#pragma link C++ class ant::TTaggerHit+;
#pragma link C++ class ant::THeaderInfo+;
#pragma link C++ class ant::TSlowControl+;
#pragma link C++ class ant::TUnpackerMessage+;

#endif // __CINT__

