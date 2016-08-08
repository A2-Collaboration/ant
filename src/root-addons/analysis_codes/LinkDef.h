#ifdef __CINT__
// some ROOTCINT defaults
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// the relevant classes
// templated classes need to have explicit type

#pragma link C++ class ant::SymmetricPi0+;
#pragma link C++ class ant::ExtractResolutions+;
#pragma link C++ class ant::CBESum_Check+;
#pragma link C++ class ant::TAPSEnergy_Check+;
#pragma link C++ class ant::DetectorPlots+;
#pragma link C++ class ant::BrowseHistogramsCanvas+;
#pragma link C++ class ant::PlotTimings+;

#pragma link C++ class ant::hstack-; // has its own Streamer implementation with cereal

#pragma link C++ class ant::histtools+;
#pragma link C++ class ant::MC::ThetaCBToyMC+;
#pragma link C++ class ant::KinfitExtract+;
#pragma link C++ class ant::Math+;
#pragma link C++ class ant::InterpolatedPulls+;
#pragma link C++ class ant::ConvergencePlot+;
#pragma link C++ class ant::DisplayClustering+;


#pragma link C++ class std::list<TDirectory*>+;

#endif // __CINT__

