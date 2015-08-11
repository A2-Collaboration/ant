#ifdef __CINT__
// some ROOTCINT defaults
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// the relevant classes
// templated classes need to have explicit type

#pragma link C++ class ant::TCalibrationEditor+;
#pragma link C++ class ant::TH2CB+;
#pragma link C++ class ant::TH2TAPS+;

#endif // __CINT__

