#ifdef __CINT__
// some ROOTCINT defaults
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// the relevant classes
// templated classes need to have explicit type

#pragma link C++ class ant::TKeyValue<std::int64_t>+;
#pragma link C++ class ant::TKeyValue<std::string>+;
#pragma link C++ class ant::TKeyValue<double>+;
#pragma link C++ class ant::TKeyValue<std::vector<double>>+;
#pragma link C++ class ant::TID+;
#pragma link C++ class ant::TDataRecord+;

#pragma link C++ class ant::TEvent-; // has its own Streamer implementation

#pragma link C++ class ant::TSlowControl+;
#pragma link C++ class ant::TUnpackerMessage+;
#pragma link C++ class ant::TCalibrationData+;

#pragma link C++ class ant::TAntHeader+;

#pragma link C++ class std::vector<TLorentzVector>+;

#endif // __CINT__

