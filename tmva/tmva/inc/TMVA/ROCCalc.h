#ifndef ROOT_TMVA_ROCCalc
#define ROOT_TMVA_ROCCalc

#include "RtypesCore.h"

class TList;
class TTree;
class TString;
class TH1;
class TH1D;
class TH2;
class TH2F;
class TSpline;

namespace TMVA {

   class MsgLogger;


   class ROCCalc {
    
   public:
      ROCCalc(TH1* mvaS, TH1* mvaB);
    
      ~ROCCalc();
    

      TH1D* GetROC();
      // return the signal eff for a given backgr. efficiency
      Double_t GetEffSForEffBof(Double_t effBref, Double_t &effSerr);
      // return the cut value 
      Double_t GetSignalReferenceCut(){return fSignalCut;} 
      // return the area under the ROC curve
      Double_t GetROCIntegral();
      // return the statistical significance as function of the mva cut value
      TH1* GetSignificance( Int_t nStot, Int_t nBtot);
      TH1* GetPurity(Int_t nStot, Int_t nBtot);
    
      void ApplySignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* any = 0 );
    
      TH1* GetMvaSpdf(){return fmvaSpdf;}
      TH1* GetMvaBpdf(){return fmvaBpdf;}
    
      //false if is found some error in mvaS or mvaB
      Bool_t GetStatus(){return fStatus;}
      void ResetStatus(){fStatus=kTRUE;}
    
   private:
      Double_t        Root(Double_t);
      Double_t        GetEffForRoot( Double_t theCut );
      Int_t           fMaxIter;  // maximum number of iterations
      Double_t        fAbsTol;   // absolute tolerance deviation
    
      Bool_t          fStatus; //false if is found some error in mvaS or mvaB 
    
      UInt_t          fNbins;
      Bool_t          fUseSplines;

      TH1*            fmvaS, *fmvaB;       // the input mva distributions
      TH1*            fmvaSpdf, *fmvaBpdf;       // the normalized (and rebinned) input mva distributions
      Float_t         fXmin, fXmax;       // min and max of the mva distribution  
      Double_t        fNevtS;             // number of signal events (used in error calculation)
      Int_t           fCutOrientation;    //+1 if larger mva value means more signal like, -1 otherwise
      TSpline*        fSplS, *fSplB;
      TSpline*        fSplmvaCumS, *fSplmvaCumB;  // spline of cumulated mva distributions
      TSpline*        fSpleffBvsS;    
      TH1*            fmvaScumul, *fmvaBcumul;
      Int_t           fnStot, fnBtot;
      TH1*            fSignificance;
      TH1*            fPurity;
      TH1D*           effBvsS;
      TH1D*           rejBvsS;
      TH1D*           inveffBvsS;

      Double_t        fSignalCut;  // MVA cut value for last demanded background rejection or signal efficiency

      mutable MsgLogger* fLogger;   //! message logger
      MsgLogger& Log() const { return *fLogger; }                       

   };
}
#endif
