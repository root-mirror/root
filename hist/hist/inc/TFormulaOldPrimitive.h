// @(#)root/hist:$Id$
// Author: Marian Ivanov, 2005

/*************************************************************************
* Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/
// ---------------------------------- TFormulaOldPrimitive.h

#ifndef ROOT_TFormulaOldPrimitive
#define ROOT_TFormulaOldPrimitive



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFormulaOldPrimitive                                                    //
//                                                                      //
// The formula primitive base class                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif
#ifndef ROOT_TBits
#include "TBits.h"
#endif
#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif

class TFormulaOld;

class TFormulaOldPrimitive : public TNamed
{
   friend class TFormulaOld;
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
      GenFuncG    fFuncG;                                 //!pointer to the TFormulaOld generic function
      GenFunc0    fFunc0;                                 //!pointer to the function
      GenFunc10   fFunc10;                                //!pointer to the function
      GenFunc110  fFunc110;                               //!pointer to the function
      GenFunc1110 fFunc1110;                              //!pointer to the function
      TFuncG      fTFuncG;                                //!pointer to the TFormulaOld generic function
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
   TFormulaOldPrimitive(const TFormulaOldPrimitive&); // Not implemented
   TFormulaOldPrimitive& operator=(const TFormulaOldPrimitive&); // Not implemented
public:
   TFormulaOldPrimitive();
   TFormulaOldPrimitive(const char *name,const char *formula, GenFunc0 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, GenFunc10 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, GenFunc110 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, GenFunc1110 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, GenFuncG fpointer,Int_t npar);
   TFormulaOldPrimitive(const char *name,const char *formula, TFunc0 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, TFunc10 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, TFunc110 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, TFunc1110 fpointer);
   TFormulaOldPrimitive(const char *name,const char *formula, TFuncG fpointer);
   static Int_t AddFormula(TFormulaOldPrimitive * formula);
   static TFormulaOldPrimitive* FindFormula(const char* name);
   static TFormulaOldPrimitive* FindFormula(const char* name, const char *args);
   static TFormulaOldPrimitive* FindFormula(const char* name, UInt_t nargs);
   Double_t Eval(Double_t* x);                   //eval primitive function
   Double_t Eval(TObject *o,  Double_t *x);      //eval member function
   Double_t Eval(Double_t *x, Double_t *param);  //eval primitive parametric function

   ClassDef(TFormulaOldPrimitive,0)  //The primitive formula
};

#endif
