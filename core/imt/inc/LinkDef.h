#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Only for the autoload, autoparse. No IO of these classes is foreseen!
// Exclude in case ROOT does not have IMT support
#ifdef R__USE_IMT
#pragma link C++ class ROOT::TThreadExecutor-;
#pragma link C++ class ROOT::Experimental::TTaskGroup-;
#endif

#endif
