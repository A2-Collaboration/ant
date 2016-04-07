#ifdef __CINT__
// some ROOTCINT defaults
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// the relevant classes
// templated classes need to have explicit type

#pragma link C++ class SymmetricPi0+;
#pragma link C++ class ExtractResolutions+;
#pragma link C++ class CBESum_Check+;
#pragma link C++ class TAPSEnergy_Check+;
#pragma link C++ class DetectorPlots+;
#pragma link C++ class BrowseHistogramsCanvas+;
#pragma link C++ class PlotTimings+;

#pragma link C++ class TestClass+;
#pragma link C++ class std::vector<TestClass*>+;
#pragma link C++ class TestClassContainer+;
#pragma link C++ class PtrTest+;

#pragma link C++ class ant::hstack-; // has its own Streamer implementation with cereal

#pragma link C++ class ant::histtools+;
#pragma link C++ class ant::MC::ThetaCBToyMC+;
#pragma link C++ class ant::KinfitExtract+;

#endif // __CINT__

