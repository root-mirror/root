// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Yair Mahalalel, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPDERS                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Multidimensional Likelihood using the "Probability density estimator      *
 *      range search" (PDERS) method suggested in                                 *
 *      T. Carli and B. Koblitz, NIM A 501, 576 (2003)                            *
 *                                                                                *
 *      The multidimensional PDFs for signal and background are modeled           *
 *      by counting the events in the "vicinity" of a test point. The volume      *
 *      that describes "vicinity" is user-defined through the option string.      *
 *      A search method based on binary-trees is used to improve the selection    *
 *      efficiency of the volume search.                                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Yair Mahalalel  <Yair.Mahalalel@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <peter.speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodPDERS
#define ROOT_TMVA_MethodPDERS

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodPDERS                                                          //
//                                                                      //
// Multidimensional Likelihood using the "Probability density           //
// estimator range search" (PDERS) method                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/MethodBase.h"
#include "TMVA/BinarySearchTree.h"
#include "TVector.h"
#include "ThreadLocalStorage.h"
#include <vector>

namespace TMVA {

   class Volume;
   class Event;

   class MethodPDERS : public MethodBase {

   public:

      MethodPDERS( const TString& jobName,
                   const TString& methodTitle,
                   DataSetInfo& theData,
                   const TString& theOption);

      MethodPDERS( DataSetInfo& theData,
                   const TString& theWeightFile);

      virtual ~MethodPDERS( void );

      virtual Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets );


      // training method
      void Train( void );

      // write weights to file
      void WriteWeightsToStream( TFile& rf ) const;
      void AddWeightsXMLTo( void* parent ) const;

      // read weights from file
      void ReadWeightsFromStream( std::istream& istr );
      void ReadWeightsFromStream( TFile& istr );
      void ReadWeightsFromXML( void* wghtnode );

      // calculate the MVA value
      Double_t GetMvaValue( Double_t* err = 0, Double_t* errUpper = 0 );

      // calculate the MVA value
      const std::vector<Float_t>& GetRegressionValues();
   public:

      // for root finder
      static Double_t IGetVolumeContentForRoot( Double_t );
      Double_t         GetVolumeContentForRoot( Double_t );

      // static pointer to this object
      static MethodPDERS* ThisPDERS( void );

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;

      Volume*      fHelpVolume; // auxiliary variable
      Int_t        fFcnCall;    // number of external function calls (RootFinder)

      // accessors
      BinarySearchTree* GetBinaryTree( void ) const { return fBinaryTree; }

      Double_t             CKernelEstimate( const Event&, std::vector<const BinarySearchTreeNode*>&, Volume& );
      void                 RKernelEstimate( const Event&, std::vector<const BinarySearchTreeNode*>&, Volume&, std::vector<Float_t> *pdfSum );

      Double_t ApplyKernelFunction( Double_t normalized_distance );
      Double_t KernelNormalization( Double_t pdf );
      Double_t GetNormalizedDistance( const TMVA::Event &base_event,
                                      const BinarySearchTreeNode &sample_event,
                                      Double_t *dim_normalization);
      Double_t NormSinc( Double_t x );
      Double_t LanczosFilter( Int_t level, Double_t x );

      // ranking of input variables
      const Ranking* CreateRanking() { return 0; }

   private:

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      // calculate the averages of the input variables needed for adaptive training
      void CalcAverages();

      // create binary search trees for signal and background
      void CreateBinarySearchTree( Types::ETreeType type );

      // get sample of training events
      void GetSample( const Event &e, std::vector<const BinarySearchTreeNode*>& events, Volume *volume);

      // option
      TString fVolumeRange;    // option volume range
      TString fKernelString;   // option kernel estimator

      enum EVolumeRangeMode {
         kUnsupported = 0,
         kMinMax,
         kRMS,
         kAdaptive,
         kUnscaled,
         kkNN
      } fVRangeMode;

      enum EKernelEstimator {
         kBox = 0,
         kSphere,
         kTeepee,
         kGauss,
         kSinc3,     // the sinc enumerators must be consecutive and in order!
         kSinc5,
         kSinc7,
         kSinc9,
         kSinc11,
         kLanczos2,
         kLanczos3,
         kLanczos5,
         kLanczos8,
         kTrim
      } fKernelEstimator;

      BinarySearchTree*  fBinaryTree;   // binary tree

      std::vector<Float_t>*   fDelta;         // size of volume
      std::vector<Float_t>*   fShift;         // volume center
      std::vector<Float_t>    fAverageRMS;    // average RMS of signal and background

      Float_t            fScaleS;        // weight for signal events
      Float_t            fScaleB;        // weight for background events
      Float_t            fDeltaFrac;     // fraction of RMS
      Double_t           fGaussSigma;    // size of Gauss in adaptive volume
      Double_t           fGaussSigmaNorm;// size of Gauss in adaptive volume (normalised to dimensions)

      Double_t           fNRegOut;       // number of output dimensions for regression

      // input for adaptive volume adjustment
      Float_t            fNEventsMin;    // minimum number of events in adaptive volume
      Float_t            fNEventsMax;    // maximum number of events in adaptive volume
      Float_t            fMaxVIterations;// maximum number of iterations to adapt volume size
      Float_t            fInitialScale;  // initial scale for adaptive volume

      Bool_t             fInitializedVolumeEle; // is volume element initialized ?

      Int_t              fkNNMin;        // min number of events in kNN tree
      Int_t              fkNNMax;        // max number of events in kNN tree

      Double_t           fMax_distance;  // maximum distance
      Bool_t             fPrinted;       // print
      Bool_t             fNormTree;      // binary-search tree is normalised

      void    SetVolumeElement ( void );

      Double_t              CRScalc           ( const Event& );
      void                  RRScalc           ( const Event&, std::vector<Float_t>* count );

      Float_t GetError         ( Float_t countS, Float_t countB,
                                 Float_t sumW2S, Float_t sumW2B ) const;

      // This is a workaround for OSx where static thread_local data members are
      // not supported. The C++ solution would indeed be the following:
      static MethodPDERS*& GetMethodPDERSThreadLocal() {TTHREAD_TLS(MethodPDERS*) fgThisPDERS(nullptr); return fgThisPDERS;};
      void UpdateThis();

      void Init( void );

      ClassDef(MethodPDERS,0); // Multi-dimensional probability density estimator range search (PDERS) method
   };

} // namespace TMVA

#endif // MethodPDERS_H
