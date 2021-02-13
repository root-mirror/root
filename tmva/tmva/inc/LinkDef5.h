#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all namespaces;

#pragma link C++ nestedclass;

#pragma link C++ namespace TMVA;

// other classes
#pragma link C++ function TMVA::annconvergencetest;                  
#pragma link C++ function TMVA::BDTControlPlots;                     
#pragma link C++ function TMVA::BDT;                                 
#pragma link C++ function TMVA::BDT_Reg;                             
#pragma link C++ function TMVA::BoostControlPlots;                   
#pragma link C++ function TMVA::compareanapp;                        
#pragma link C++ function TMVA::correlationscatters;                 
#pragma link C++ function TMVA::correlationscattersMultiClass;       
#pragma link C++ function TMVA::correlations;                        
#pragma link C++ function TMVA::correlationsMultiClass;              
#pragma link C++ function TMVA::CorrGui;                             
#pragma link C++ function TMVA::CorrGuiMultiClass;                   
#pragma link C++ function TMVA::deviations;                          
#pragma link C++ function TMVA::efficiencies;                        
#pragma link C++ function TMVA::likelihoodrefs;                      
#pragma link C++ function TMVA::MovieMaker;                          
#pragma link C++ defined_in "TMVA/mvaeffs.h";                             
#pragma link C++ class    TMVA::MethodInfo+;                                
#pragma link C++ function TMVA::mvas;                                
#pragma link C++ function TMVA::mvasMulticlass;                      
#pragma link C++ function TMVA::mvaweights;                          
#pragma link C++ function TMVA::network;                             
#pragma link C++ function TMVA::paracoor;                            
//#pragma link C++ function TMVA::PlotFoams;                           
#pragma link C++ defined_in "TMVA/PlotFoams.h";                           
#pragma link C++ function TMVA::probas;                              
#pragma link C++ function TMVA::regression_averagedevs;              
#pragma link C++ function TMVA::rulevisCorr;                         
#pragma link C++ function TMVA::rulevis;                             
#pragma link C++ function TMVA::rulevisHists;                        
#pragma link C++ function TMVA::TMVAGui;                             
#pragma link C++ function TMVA::TMVAMultiClassGui;                   
#pragma link C++ function TMVA::TMVARegGui;                          
#pragma link C++ function TMVA::variables;                           
#pragma link C++ function TMVA::variablesMultiClass;   

#pragma link C++ class TMVA::StatDialogBDT+;
#pragma link C++ class TMVA::StatDialogBDTReg+;
#pragma link C++ class TMVA::StatDialogMVAEffs+;


#endif
