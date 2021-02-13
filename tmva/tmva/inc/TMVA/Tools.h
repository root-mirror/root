// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Tools                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Global auxiliary applications and data treatment routines                 *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <peter.speckmayer@cern.ch> - CERN, Switzerland           *
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

#ifndef ROOT_TMVA_Tools
#define ROOT_TMVA_Tools

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Tools (namespace)                                                    //
//                                                                      //
// Global auxiliary applications and data treatment routines            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#if __cplusplus > 199711L
#include <atomic>
#endif

#include "TXMLEngine.h"

#include "TMatrixDSymfwd.h"

#include "TMatrixDfwd.h"

#include "TVectorDfwd.h"

#include "TMVA/Types.h"

#include "TMVA/VariableTransformBase.h"

#include "TString.h"

#include "TMVA/MsgLogger.h"

class TList;
class TTree;
class TH1;
class TH2;
class TH2F;
class TSpline;
class TXMLEngine;

namespace TMVA {

   class Event;
   class PDF;
   class MsgLogger;

   class Tools {

   private:

      Tools();

   public:

      // destructor
      ~Tools();

      // accessor to single instance
      static Tools& Instance();
      static void   DestroyInstance();


      template <typename T> Double_t Mean(Long64_t n, const T *a, const Double_t *w=0);
      template <typename Iterator, typename WeightIterator> Double_t Mean ( Iterator first, Iterator last, WeightIterator w);

      template <typename T> Double_t RMS(Long64_t n, const T *a, const Double_t *w=0);
      template <typename Iterator, typename WeightIterator> Double_t RMS(Iterator first, Iterator last, WeightIterator w);


      // simple statistics operations on tree entries
      void  ComputeStat( const std::vector<TMVA::Event*>&,
                         std::vector<Float_t>*,
                         Double_t&, Double_t&, Double_t&,
                         Double_t&, Double_t&, Double_t&, Int_t signalClass,
                         Bool_t norm = kFALSE );

      // compute variance from sums
      inline Double_t ComputeVariance( Double_t sumx2, Double_t sumx, Int_t nx );

      // creates histograms normalized to one
      TH1* projNormTH1F( TTree* theTree, const TString& theVarName,
                         const TString& name, Int_t nbins,
                         Double_t xmin, Double_t xmax, const TString& cut );

      // normalize histogram by its integral
      Double_t NormHist( TH1* theHist, Double_t norm = 1.0 );

      // parser for TString phrase with items separated by a character
      TList* ParseFormatLine( TString theString, const char * sep = ":" );

      // parse option string for ANN methods
      std::vector<Int_t>* ParseANNOptionString( TString theOptions, Int_t nvar,
                                                std::vector<Int_t>* nodes );

      // returns the square-root of a symmetric matrix: symMat = sqrtMat*sqrtMat
      TMatrixD* GetSQRootMatrix( TMatrixDSym* symMat );

      // returns the covariance matrix of of the different classes (and the sum)
      // given the event sample
      std::vector<TMatrixDSym*>* CalcCovarianceMatrices( const std::vector<Event*>& events, Int_t maxCls, VariableTransformBase* transformBase=0 );
      std::vector<TMatrixDSym*>* CalcCovarianceMatrices( const std::vector<const Event*>& events, Int_t maxCls, VariableTransformBase* transformBase=0 );


      // turns covariance into correlation matrix
      const TMatrixD* GetCorrelationMatrix( const TMatrixD* covMat );

      // check spline quality by comparison with initial histogram
      Bool_t CheckSplines( const TH1*, const TSpline* );

      // normalization of variable output
      Double_t NormVariable( Double_t x, Double_t xmin, Double_t xmax );

      // return separation of two histograms
      Double_t GetSeparation( TH1* S, TH1* B ) const;
      Double_t GetSeparation( const PDF& pdfS, const PDF& pdfB ) const;

      // vector rescaling
      std::vector<Double_t> MVADiff( std::vector<Double_t>&, std::vector<Double_t>& );
      void Scale( std::vector<Double_t>&, Double_t );
      void Scale( std::vector<Float_t>&,  Float_t  );

      // re-arrange a vector of arrays (vectors) in a way such that the first array
      // is ordered, and the other arrays reshuffled accordingly
      void UsefulSortDescending( std::vector< std::vector<Double_t> >&, std::vector<TString>* vs = 0 );
      void UsefulSortAscending ( std::vector< std::vector<Double_t> >&, std::vector<TString>* vs = 0 );

      void UsefulSortDescending( std::vector<Double_t>& );
      void UsefulSortAscending ( std::vector<Double_t>& );

      Int_t GetIndexMaxElement ( std::vector<Double_t>& );
      Int_t GetIndexMinElement ( std::vector<Double_t>& );

      // check if input string contains regular expression
      Bool_t  ContainsRegularExpression( const TString& s );
      TString ReplaceRegularExpressions( const TString& s, const TString& replace = "+" );

      // routines for formatted output -----------------
      void FormattedOutput( const std::vector<Double_t>&, const std::vector<TString>&,
                            const TString titleVars, const TString titleValues, MsgLogger& logger,
                            TString format = "%+1.3f" );
      void FormattedOutput( const TMatrixD&, const std::vector<TString>&, MsgLogger& logger );
      void FormattedOutput( const TMatrixD&, const std::vector<TString>& vert, const std::vector<TString>& horiz,
                            MsgLogger& logger );

      void WriteFloatArbitraryPrecision( Float_t  val, std::ostream& os );
      void ReadFloatArbitraryPrecision ( Float_t& val, std::istream& is );

      // for histogramming
      TString GetXTitleWithUnit( const TString& title, const TString& unit );
      TString GetYTitleWithUnit( const TH1& h, const TString& unit, Bool_t normalised );

      // Mutual Information method for non-linear correlations estimates in 2D histogram
      // Author: Moritz Backes, Geneva (2009)
      Double_t GetMutualInformation( const TH2F& );

      // Correlation Ratio method for non-linear correlations estimates in 2D histogram
      // Author: Moritz Backes, Geneva (2009)
      Double_t GetCorrelationRatio( const TH2F& );
      TH2F*    TransposeHist      ( const TH2F& );

      // check if "silent" or "verbose" option in configuration string
      Bool_t CheckForSilentOption ( const TString& ) const;
      Bool_t CheckForVerboseOption( const TString& ) const;

      // color information
      const TString& Color( const TString& );

      // print welcome message (to be called from, eg, .TMVAlogon)
      enum EWelcomeMessage { kStandardWelcomeMsg = 1,
                             kIsometricWelcomeMsg,
                             kBlockWelcomeMsg,
                             kLeanWelcomeMsg,
                             kLogoWelcomeMsg,
                             kSmall1WelcomeMsg,
                             kSmall2WelcomeMsg,
                             kOriginalWelcomeMsgColor,
                             kOriginalWelcomeMsgBW };

      // print TMVA citation (to be called from, eg, .TMVAlogon)
      enum ECitation { kPlainText = 1,
                       kBibTeX,
                       kLaTeX,
                       kHtmlLink };

      void TMVAWelcomeMessage();
      void TMVAWelcomeMessage( MsgLogger& logger, EWelcomeMessage m = kStandardWelcomeMsg );
      void TMVAVersionMessage( MsgLogger& logger );
      void ROOTVersionMessage( MsgLogger& logger );

      void TMVACitation( MsgLogger& logger, ECitation citType = kPlainText );

      // string tools

      std::vector<TString> SplitString( const TString& theOpt, const char separator ) const;

      // variables
      const TString fRegexp;
      mutable MsgLogger*    fLogger;
      MsgLogger& Log() const { return *fLogger; }
#if __cplusplus > 199711L
      static std::atomic<Tools*> fgTools;
#else
      static Tools* fgTools;
#endif

      // xml tools

      TString     StringFromInt      ( Long_t i   );
      TString     StringFromDouble   ( Double_t d );
      void        WriteTMatrixDToXML ( void* node, const char* name, TMatrixD* mat );
      void        WriteTVectorDToXML ( void* node, const char* name, TVectorD* vec );
      void        ReadTMatrixDFromXML( void* node, const char* name, TMatrixD* mat );
      void        ReadTVectorDFromXML( void* node, const char* name, TVectorD* vec );
      Bool_t      HistoHasEquidistantBins(const TH1& h);

      Bool_t      HasAttr     ( void* node, const char* attrname );
      template<typename T>
         inline void ReadAttr    ( void* node, const char* , T& value );
      void        ReadAttr    ( void* node, const char* attrname, TString& value );
      void ReadAttr(void *node, const char *, float &value);
      void ReadAttr(void *node, const char *, int &value);
      void ReadAttr(void *node, const char *, short &value);

      template<typename T>
         void        AddAttr     ( void* node, const char* , const T& value, Int_t precision = 16 );
      void        AddAttr     ( void* node, const char* attrname, const char* value );
      void*       AddChild    ( void* parent, const char* childname, const char* content = 0, bool isRootNode = false );
      Bool_t      AddRawLine  ( void* node, const char * raw );
      Bool_t      AddComment  ( void* node, const char* comment );

      void*       GetParent( void* child);
      void*       GetChild    ( void* parent, const char* childname=0 );
      void*       GetNextChild( void* prevchild, const char* childname=0 );
      const char* GetContent  ( void* node );
      const char* GetName     ( void* node );

      TXMLEngine& xmlengine() { return *fXMLEngine; }
      int xmlenginebuffersize() { return fXMLBufferSize;}
      void SetXMLEngineBufferSize(int buffer) { fXMLBufferSize = buffer; }
      TXMLEngine* fXMLEngine;

      TH1*       GetCumulativeDist( TH1* h);

   private:

      int fXMLBufferSize = 10000000; 
      // utilities for correlation ratio
      Double_t GetYMean_binX( const TH2& , Int_t bin_x );

   }; // Common tools

   Tools& gTools(); // global accessor

   //
   // Adapts a TRandom random number generator to the interface of the ones in the
   // standard library (STL) so that TRandom derived generators can be used with
   // STL algorithms such as `std::shuffle`.
   //
   // Example:
   // ```
   // std::vector<double> v {0, 1, 2, 3, 4, 5};
   // TRandom3StdEngine rng(seed);
   // std::shuffle(v.begin(), v.end(), rng);
   // ```
   //
   // Or at a lower level:
   // ```
   // std::vector<double> v {0, 1, 2, 3, 4, 5};
   // RandomGenerator<TRandom3> rng(seed);
   // std::shuffle(v.begin(), v.end(), rng);
   // ```
   //
   template <typename TRandomLike, typename UIntType = UInt_t, UIntType max_val = kMaxUInt>
   class RandomGenerator {
   public:
      using result_type = UIntType;

      RandomGenerator(UIntType s = 0) { fRandom.SetSeed(s); }

      static constexpr UIntType min() { return 0; }
      static constexpr UIntType max() { return max_val; }

      void seed(UIntType s = 0) { fRandom.SetSeed(s); }

      UIntType operator()() { return fRandom.Integer(max()); }

      void discard(unsigned long long z)
      {
         double r;
         for (unsigned long long i = 0; i < z; ++i)
            r = fRandom.Rndm();
         (void) r; /* avoid unused variable warning */
      }

   private:
      TRandomLike fRandom; // random generator
   };

} // namespace TMVA

////////////////////////////////////////////////////////////////////////////////
/// read attribute from xml

template<typename T> void TMVA::Tools::ReadAttr( void* node, const char* attrname, T& value )
{
   // read attribute from xml
   const char *val = xmlengine().GetAttr(node, attrname);
   if (val == 0) {
      const char *nodename = xmlengine().GetNodeName(node);
      Log() << kFATAL << "Trying to read non-existing attribute '" << attrname << "' from xml node '" << nodename << "'"
            << Endl;
   }
   std::stringstream s(val);
   // coverity[tainted_data_argument]
   s >> value;
}

////////////////////////////////////////////////////////////////////////////////
/// add attribute to xml

template<typename T>
void TMVA::Tools::AddAttr( void* node, const char* attrname, const T& value, Int_t precision )
{
   std::stringstream s;
   s.precision( precision );
   s << std::scientific << value;
   AddAttr( node, attrname, s.str().c_str() );
}

////////////////////////////////////////////////////////////////////////////////
/// compute variance from given sums

inline Double_t TMVA::Tools::ComputeVariance( Double_t sumx2, Double_t sumx, Int_t nx )
{
   if (nx<2) return 0;
   return (sumx2 - ((sumx*sumx)/static_cast<Double_t>(nx)))/static_cast<Double_t>(nx-1);
}

#endif
