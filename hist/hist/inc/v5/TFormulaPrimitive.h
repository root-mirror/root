// @(#)root/hist:$Id$
// Author: Marian Ivanov, 2005

/*************************************************************************
* Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/
// ---------------------------------- TFormulaPrimitive.h

#ifndef ROOT_v5_TFormulaPrimitive
#define ROOT_v5_TFormulaPrimitive



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFormulaPrimitive                                                    //
//                                                                      //
// The formula primitive base class                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TBits.h"

class TObjArray;

namespace ROOT {

   namespace v5 {

class TFormula;

class TFormulaPrimitive : public TNamed
{
   friend class ROOT::v5::TFormula;
public:
   typedef Double_t (*GenFuncG)(const Double_t*,const Double_t*);
   typedef Double_t (*GenFunc0)();
   typedef Double_t (*GenFunc10)(Double_t);
   typedef Double_t (*GenFunc110)(Double_t,Double_t);
   typedef Double_t (*GenFunc1110)(Double_t,Double_t, Double_t);
   typedef Double_t (TObject::*TFuncG)(const Double_t*,const Double_t*) const;
   typedef Double_t (TObject::*TFunc0)() const;
   typedef Double_t (TObject::*TFunc10)(Double_t) const;
   typedef Double_t (TObject::*TFunc110)(Double_t,Double_t) const;
   typedef Double_t (TObject::*TFunc1110)(Double_t,Double_t,Double_t) const;
protected:
   static TObjArray * fgListOfFunction;                   //!list of global primitive formulas
   static Int_t       BuildBasicFormulas();               //build list of basic formulas
   union {
      GenFuncG    fFuncG;                                 //!pointer to the TFormula generic function
      GenFunc0    fFunc0;                                 //!pointer to the function
      GenFunc10   fFunc10;                                //!pointer to the function
      GenFunc110  fFunc110;                               //!pointer to the function
      GenFunc1110 fFunc1110;                              //!pointer to the function
      TFuncG      fTFuncG;                                //!pointer to the TFormula generic function
      TFunc0      fTFunc0;                                //! pointer to member function
      TFunc10     fTFunc10;                               //! pointer to member function
      TFunc110    fTFunc110;                              //! pointer to member function
      TFunc1110   fTFunc1110;                             //! pointer to member function
   };
   Int_t      fType;                                      //type of the function
   Int_t      fNArguments;                                //number of arguments
   Int_t      fNParameters;                               //number of parameters
   Bool_t     fIsStatic;                                  // indication if the function is static
private:
   TFormulaPrimitive(const TFormulaPrimitive&); // Not implemented
   TFormulaPrimitive& operator=(const TFormulaPrimitive&); // Not implemented
public:
   TFormulaPrimitive();
   TFormulaPrimitive(const char *name,const char *formula, GenFunc0 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, GenFunc10 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, GenFunc110 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, GenFunc1110 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, GenFuncG fpointer,Int_t npar);
   TFormulaPrimitive(const char *name,const char *formula, TFunc0 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, TFunc10 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, TFunc110 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, TFunc1110 fpointer);
   TFormulaPrimitive(const char *name,const char *formula, TFuncG fpointer);
   static Int_t AddFormula(TFormulaPrimitive * formula);
   static TFormulaPrimitive* FindFormula(const char* name);
   static TFormulaPrimitive* FindFormula(const char* name, const char *args);
   static TFormulaPrimitive* FindFormula(const char* name, UInt_t nargs);
   Double_t Eval(Double_t* x);                   //eval primitive function
   Double_t Eval(TObject *o,  Double_t *x);      //eval member function
   Double_t Eval(Double_t *x, Double_t *param);  //eval primitive parametric function

   ClassDef(ROOT::v5::TFormulaPrimitive,0)  //The primitive formula
};

   } // end namespace v5

} // end namespace ROOT


#endif
