// @(#)root/tmva $Id$
// Author: S. Jadach, Tancredi Carli, Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoam                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Class for PDEFoam object                                                  *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      S. Jadach        - Institute of Nuclear Physics, Cracow, Poland           *
 *      Tancredi Carli   - CERN, Switzerland                                      *
 *      Dominik Dannheim - CERN, Switzerland                                      *
 *      Alexander Voigt  - TU Dresden, Germany                                    *
 *                                                                                *
 * Copyright (c) 2008, 2010:                                                      *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_PDEFoam
#define ROOT_TMVA_PDEFoam

#include <iosfwd>
#include <cassert>
#include <vector>
#include <map>

#include "TH2D.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TVectorT.h"
#include "TString.h"
#include "TMVA/VariableInfo.h"
#include "TMVA/Timer.h"
#include "TObject.h"
#include "TRandom3.h"

namespace TMVA {
   class PDEFoamCell;
   class PDEFoamVect;
   class PDEFoamDensityBase;
   class PDEFoamKernelBase;
   class PDEFoam;

   // separation types
   enum EDTSeparation { kFoam, kGiniIndex, kMisClassificationError,
                        kCrossEntropy, kGiniIndexWithLaplace, kSdivSqrtSplusB };

   // foam types
   enum EFoamType { kSeparate, kDiscr, kMonoTarget, kMultiTarget, kMultiClass };

   // enum type for possible foam cell values
   // kValue         : cell value who's rms is minimized
   // kValueError    : error on kValue
   // kValueDensity  : kValue / cell volume
   // kMeanValue     : mean sampling value (saved in fIntegral)
   // kRms           : rms of sampling distribution (saved in fDriver)
   // kRmsOvMean     : rms/mean of sampling distribution (saved in
   //                  fDriver and fIntegral)
   // kCellVolume    : volume of cell
   enum ECellValue { kValue, kValueError, kValueDensity, kMeanValue,
                     kRms, kRmsOvMean, kCellVolume };
}

#include "TMVA/PDEFoamDensityBase.h"
#include "TMVA/PDEFoamVect.h"
#include "TMVA/PDEFoamCell.h"

namespace TMVA {

   class PDEFoam : public TObject {
   protected:
      // COMPONENTS //
      //-------------- Input parameters
      TString fName;             // Name of a given instance of the FOAM class
      Int_t   fDim;              // Dimension of the integration/simulation space
      Int_t   fNCells;           // Maximum number of cells
      //-------------------
      Int_t   fNBin;             // No. of bins in the edge histogram for cell MC exploration
      Int_t   fNSampl;           // No. of MC events, when dividing (exploring) cell
      Int_t   fEvPerBin;         // Maximum number of effective (wt=1) events per bin
      //-------------------  MULTI-BRANCHING ---------------------
      Int_t  *fMaskDiv;          //! [fDim] Dynamic Mask for cell division
      Int_t  *fInhiDiv;          //! [fDim] Flags for inhibiting cell division
      //-------------------  GEOMETRY ----------------------------
      Int_t   fNoAct;            // Number of active cells
      Int_t   fLastCe;           // Index of the last cell
      PDEFoamCell **fCells;      // [fNCells] Array of ALL cells
      //------------------ M.C. generation----------------------------
      TObjArray *fHistEdg;       // Histograms of wt, one for each cell edge
      Double_t *fRvec;           // [fDim] random number vector from r.n. generator fDim+1 maximum elements
      //----------- Procedures
      TRandom3        *fPseRan;  // Pointer to user-defined generator of pseudorandom numbers
      //----------  working space for CELL exploration -------------
      Double_t *fAlpha;          // [fDim] Internal parameters of the hyperrectangle
      // ---------  PDE-Foam specific variables
      EFoamType fFoamType;     // BACKWARDS COMPATIBILITY: type of foam
      Double_t *fXmin;         // [fDim] minimum for variable transform
      Double_t *fXmax;         // [fDim] maximum for variable transform
      UInt_t fNElements;       // BACKWARDS COMPATIBILITY: number of variables in every cell
      UInt_t fNmin;            // minimal number of events in cell to split cell
      UInt_t fMaxDepth;        // maximum depth of cell tree
      Float_t fVolFrac;        // BACKWARDS COMPATIBILITY: volume fraction (with respect to total phase space
      Bool_t fFillFoamWithOrigWeights; // BACKWARDS COMPATIBILITY: fill the foam with boost or orig. weights
      EDTSeparation fDTSeparation; // BACKWARDS COMPATIBILITY: split cells according to decision tree logic
      Bool_t fPeekMax;         // BACKWARDS COMPATIBILITY: peek up cell with max. driver integral for split
      PDEFoamDensityBase *fDistr;  //! distribution of training events
      Timer *fTimer;           //! timer for graphical output
      TObjArray *fVariableNames;// collection of all variable names
      mutable MsgLogger* fLogger;                     //! message logger

      /////////////////////////////////////////////////////////////////
      //                            METHODS                          //
      /////////////////////////////////////////////////////////////////

   protected:
      // ---------- TMVA console output

      void OutputGrow(Bool_t finished = false ); // nice TMVA console output

      // ---------- Foam build-up functions

      // Internal foam initialization functions
      void InitCells();                   // Initialisation of all foam cells
      Int_t CellFill(Int_t, PDEFoamCell*);// Allocates new empty cell and return its index
      virtual void Explore(PDEFoamCell *Cell); // Exploration of the new cell, determine <wt>, wtMax etc.
      void Varedu(Double_t [], Int_t&, Double_t&,Double_t&); // Determines the best edge, variance reduction
      void MakeAlpha();             // Provides random point inside hyperrectangle
      void Grow();                  // build up foam
      Long_t PeekMax();             // peek cell with max. driver integral
      Int_t  Divide(PDEFoamCell *); // Divide iCell into two daughters; iCell retained, taged as inactive
      Double_t Eval(Double_t *xRand, Double_t &event_density); // evaluate distribution on point 'xRand'

      // ---------- Cell value access functions

      // low level functions to access a certain cell value
      Double_t GetCellElement(const PDEFoamCell *cell, UInt_t i) const;  // get Element 'i' in cell 'cell'
      void SetCellElement(PDEFoamCell *cell, UInt_t i, Double_t value); // set Element 'i' in cell 'cell' to value 'value'

      // specific function used during evaluation; determines, whether a cell value is undefined
      virtual Bool_t  CellValueIsUndefined( PDEFoamCell* );

      // finds cell according to given event variables
      PDEFoamCell* FindCell(const std::vector<Float_t>&) const ;
      std::vector<TMVA::PDEFoamCell*> FindCells(const std::vector<Float_t>&) const;
      std::vector<TMVA::PDEFoamCell*> FindCells(const std::map<Int_t,Float_t>&) const;
      void FindCells(const std::map<Int_t, Float_t>&, PDEFoamCell*, std::vector<PDEFoamCell*> &) const;

      // get internal density
      PDEFoamDensityBase* GetDistr() const { assert(fDistr); return fDistr; }

      // Square function (fastest implementation)
      template<typename T> T Sqr(T x) const { return x*x; }

      PDEFoam(const PDEFoam&);    // Copy Constructor  NOT USED

      // ---------- Public functions ----------------------------------
   public:
      PDEFoam();                  // Default constructor (used only by ROOT streamer)
      PDEFoam(const TString&);    // Principal user-defined constructor
      virtual ~PDEFoam();         // Default destructor

      // ---------- Foam creation functions

      void Initialize() {}        // initialize the PDEFoam
      void FillBinarySearchTree( const Event* ev ); // fill event into BST
      void Create();              // build-up foam

      // function to fill created cell with given value
      virtual void FillFoamCells(const Event* ev, Float_t wt);

      // remove all cell elements
      void ResetCellElements();

      // function to call after foam is grown
      virtual void Finalize() {}

      // ---------- Getters and Setters

      void SetDim(Int_t kDim); // Sets dimension of cubical space
      void SetnCells(Long_t nCells){fNCells =nCells;}  // Sets maximum number of cells
      void SetnSampl(Long_t nSampl){fNSampl =nSampl;}  // Sets no of MC events in cell exploration
      void SetnBin(Int_t nBin){fNBin = nBin;}          // Sets no of bins in histograms in cell exploration
      void SetEvPerBin(Int_t EvPerBin){fEvPerBin =EvPerBin;} // Sets max. no. of effective events per bin
      void SetInhiDiv(Int_t, Int_t ); // Set inhibition of cell division along certain edge
      void SetDensity(PDEFoamDensityBase *dens) { fDistr = dens; }

      // coverity[ -tainted_data_return ]
      Int_t    GetTotDim()    const {return fDim;  } // Get total dimension
      TString  GetFoamName()  const {return fName; } // Get name of foam
      UInt_t   GetNActiveCells()   const {return fNoAct;} // returns number of active cells
      UInt_t   GetNInActiveCells() const {return GetNCells()-GetNActiveCells();} // returns number of not active cells
      UInt_t   GetNCells()         const {return fNCells;}   // returns number of cells
      PDEFoamCell* GetRootCell()   const {return fCells[0];} // get pointer to root cell

      // Getters and Setters for user cut options
      void     SetNmin(UInt_t val)     { fNmin=val;      }
      UInt_t   GetNmin()               { return fNmin;   }
      void     SetMaxDepth(UInt_t maxdepth) { fMaxDepth = maxdepth; }
      UInt_t   GetMaxDepth() const { return fMaxDepth; }

      // Getters and Setters for foam boundaries
      void SetXmin(Int_t idim, Double_t wmin);
      void SetXmax(Int_t idim, Double_t wmax);
      Double_t GetXmin(Int_t idim) const {return fXmin[idim];}
      Double_t GetXmax(Int_t idim) const {return fXmax[idim];}

      // Getters and Setters for variable names
      void AddVariableName(const char *s) { AddVariableName(new TObjString(s)); }
      void AddVariableName(TObjString *s) { fVariableNames->Add(s); }
      TObjString* GetVariableName(Int_t idx) {return dynamic_cast<TObjString*>(fVariableNames->At(idx));}

      // Delete the fDistr object, which contains the binary search
      // tree
      void DeleteBinarySearchTree();

      // ---------- Transformation functions for event variables into foam boundaries
      // reason: foam always has boundaries [0, 1]

      Float_t VarTransform(Int_t idim, Float_t x) const; // transform [xmin, xmax] --> [0, 1]
      std::vector<Float_t> VarTransform(const std::vector<Float_t> &invec) const;
      Float_t VarTransformInvers(Int_t idim, Float_t x) const; // transform [0, 1] --> [xmin, xmax]
      std::vector<Float_t> VarTransformInvers(const std::vector<Float_t> &invec) const;

      // ---------- Debug functions

      void     CheckAll(Int_t);  // Checks correctness of the entire data structure in the FOAM object
      void     PrintCell(Long_t iCell=0); // Print content of cell
      void     PrintCells();     // Prints content of all cells

      // Message logger
      MsgLogger& Log() const { return *fLogger; }

      // ---------- Foam projection methods

      // project foam to two-dimensional histogram
      virtual TH2D* Project2(Int_t idim1, Int_t idim2, ECellValue cell_value=kValue,
                             PDEFoamKernelBase *kernel=NULL, UInt_t nbin=50);

      // Project one-dimensional foam to a 1-dim histogram
      TH1D* Draw1Dim(ECellValue cell_value, Int_t nbin, PDEFoamKernelBase *kernel=NULL);

      // Generates C++ code (root macro) for drawing foam with boxes (only 2-dim!)
      void RootPlot2dim( const TString& filename, TString opt,
                         Bool_t createCanvas = kTRUE, Bool_t colors = kTRUE );

      // ---------- Foam evaluation functions

      // get cell value for a given event
      virtual Float_t GetCellValue( const std::vector<Float_t>& xvec, ECellValue cv, PDEFoamKernelBase*  );

      // get cell values for a given (incomplete) event vector
      virtual std::vector<Float_t> GetCellValue( const std::map<Int_t,Float_t>& xvec, ECellValue cv );

      // get cell value stored in a foam cell
      virtual Float_t GetCellValue( const PDEFoamCell* cell, ECellValue cv );

      // ---------- friend classes
      friend class PDEFoamKernelBase;
      friend class PDEFoamKernelTrivial;
      friend class PDEFoamKernelLinN;
      friend class PDEFoamKernelGauss;

      // ---------- ROOT class definition
      ClassDef(PDEFoam,7) // Tree of PDEFoamCells
         }; // end of PDEFoam

}  // namespace TMVA

// ---------- Inline functions

//_____________________________________________________________________
inline Float_t TMVA::PDEFoam::VarTransform(Int_t idim, Float_t x) const
{
   // transform variable x from [xmin, xmax] --> [0, 1]
   return (x-fXmin[idim])/(fXmax[idim]-fXmin[idim]);
}

//_____________________________________________________________________
inline std::vector<Float_t> TMVA::PDEFoam::VarTransform(const std::vector<Float_t> &invec) const
{
   // transform vector invec from [xmin, xmax] --> [0, 1]
   std::vector<Float_t> outvec;
   for(UInt_t i=0; i<invec.size(); i++)
      outvec.push_back(VarTransform(i, invec.at(i)));
   return outvec;
}

//_____________________________________________________________________
inline Float_t TMVA::PDEFoam::VarTransformInvers(Int_t idim, Float_t x) const
{
   // transform variable x from [0, 1] --> [xmin, xmax]
   return x*(fXmax[idim]-fXmin[idim]) + fXmin[idim];
}

//_____________________________________________________________________
inline std::vector<Float_t> TMVA::PDEFoam::VarTransformInvers(const std::vector<Float_t> &invec) const
{
   // transform vector invec from [0, 1] --> [xmin, xmax]
   std::vector<Float_t> outvec;
   for(UInt_t i=0; i<invec.size(); i++)
      outvec.push_back(VarTransformInvers(i, invec.at(i)));
   return outvec;
}

#endif
