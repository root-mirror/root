/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooExtendedBinding.h" 
#include "RooAbsPdf.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooExtendedBinding); 

 RooExtendedBinding::RooExtendedBinding(const char *name, const char *title, RooAbsPdf& _pdf) :
   RooAbsReal(name,title), 
   pdf("pdf","pdf",this,_pdf)
 { 
 } 


 RooExtendedBinding::RooExtendedBinding(const RooExtendedBinding& other, const char* name) :  
   RooAbsReal(other,name), 
   pdf("pdf",this,other.pdf)
 { 
 } 



 Double_t RooExtendedBinding::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
   return ((RooAbsPdf &)pdf.arg()).expectedEvents(nullptr);
 } 



