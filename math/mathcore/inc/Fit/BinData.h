// @(#)root/mathcore:$Id$
// Author: L. Moneta Wed Aug 30 11:15:23 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class BinData

#ifndef ROOT_Fit_BinData
#define ROOT_Fit_BinData

#include "Fit/FitData.h"
#include "Math/Error.h"
#include <cmath>
#include <vector>

namespace ROOT {

   namespace Fit {



//___________________________________________________________________________________
/**
   Class describing the binned data sets :
   vectors of  x coordinates, y values and optionally error on y values and error on coordinates
   The dimension of the coordinate is free
   There are 4 different options:
   - only coordinates and values  (for binned likelihood fits)  : kNoError
   - coordinate, values and error on  values (for normal least square fits)  : kValueError
   - coordinate, values, error on values and coordinates (for effective least square fits) : kCoordError
   - corrdinate, values, error on coordinates and asymmettric error on valyes : kAsymError

   In addition there is the option to construct Bindata copying the data in (using the DataVector class)
   or using pointer to external data (DataWrapper) class.
   In general is found to be more efficient to copy the data.
   In case of really large data sets for limiting memory consumption then the other option can be used
   Specialized constructor exists for data up to 3 dimensions.

   When the data are copying in the number of points can be set later (or re-set) using Initialize and
   the data are inserted one by one using the Add method.
   It is mandatory to set the size before using the Add method.

   @ingroup  FitData
*/


class BinData  : public FitData  {

public :

   enum ErrorType { kNoError, kValueError, kCoordError, kAsymError };

   /**
      constructor from dimension of point  and max number of points (to pre-allocate vector)
      Give a zero value and then use Initialize later one if the size is not known
   */

   explicit BinData(unsigned int maxpoints = 0, unsigned int dim = 1,
                    ErrorType err = kValueError);


   /**
      constructor from option and default range
   */
   explicit BinData (const DataOptions & opt, unsigned int maxpoints = 0,
                     unsigned int dim = 1, ErrorType err = kValueError);

   /**
      constructor from options and range
      efault is 1D and value errors
   */
   BinData (const DataOptions & opt, const DataRange & range,
            unsigned int maxpoints = 0, unsigned int dim = 1, ErrorType err = kValueError );

   /** constructurs using external data */

   /**
      constructor from external data for 1D with errors on  coordinate and value
   */
   BinData(unsigned int n, const double * dataX, const double * val,
           const double * ex , const double * eval );

   /**
      constructor from external data for 2D with errors on  coordinate and value
   */
   BinData(unsigned int n, const double * dataX, const double * dataY,
           const double * val, const double * ex , const double * ey,
           const double * eval  );

   /**
      constructor from external data for 3D with errors on  coordinate and value
   */
   BinData(unsigned int n, const double * dataX, const double * dataY,
           const double * dataZ, const double * val, const double * ex ,
           const double * ey , const double * ez , const double * eval   );

   /**
      destructor
   */
   virtual ~BinData();

   /**
      copy constructors
   */
   BinData(const BinData & rhs);

   BinData & operator= ( const BinData & rhs );


   /**
      preallocate a data set with given size ,  dimension and error type (to get the full point size)
      If the data set already exists and it is having the compatible point size space for the new points
      is created in the data sets, while if not compatible the old data are erased and new space of
      new size is allocated.
      (i.e if exists initialize is equivalent to a resize( NPoints() + maxpoints)
   */

   void Append( unsigned int newPoints, unsigned int dim = 1, ErrorType err = kValueError );

   void Initialize( unsigned int newPoints, unsigned int dim = 1, ErrorType err = kValueError );

   /**
      flag to control if data provides error on the coordinates
   */
   bool HaveCoordErrors() const {
      assert (  fErrorType == kNoError ||
                fErrorType == kValueError ||
                fErrorType == kCoordError ||
                fErrorType == kAsymError );

      return fErrorType == kCoordError;
   }

   /**
      flag to control if data provides asymmetric errors on the value
   */
   bool HaveAsymErrors() const {
      assert (  fErrorType == kNoError ||
                fErrorType == kValueError ||
                fErrorType == kCoordError ||
                fErrorType == kAsymError );

      return fErrorType == kAsymError;
   }


   /**
      apply a Log transformation of the data values
      can be used for example when fitting an exponential or gaussian
      Transform the data in place need to copy if want to preserve original data
      The data sets must not contain negative values. IN case it does,
      an empty data set is returned
   */
   BinData & LogTransform();


   /**
      add one dim data with only coordinate and values
   */
   void Add( double x, double y );

   /**
      add one dim data with no error in the coordinate (x)
      in this case store the inverse of the error in the value (y)
   */
   void Add( double x, double y, double ey );

   /**
      add one dim data with  error in the coordinate (x)
      in this case store the value (y)  error and not the inverse
   */
   void Add( double x, double y, double ex, double ey );

   /**
      add one dim data with  error in the coordinate (x) and asymmetric errors in the value (y)
      in this case store the y errors and not the inverse
   */
   void Add( double x, double y, double ex, double eyl, double eyh );

   /**
      add multi-dim coordinate data with only value
   */
   void Add( const double* x, double val );

   /**
      add multi-dim coordinate data with only error in value
   */
   void Add( const double* x, double val, double eval );

   /**
      add multi-dim coordinate data with both error in coordinates and value
   */
   void Add( const double* x, double val, const double* ex, double eval );

   /**
      add multi-dim coordinate data with both error in coordinates and value
   */
   void Add( const double* x, double val, const double* ex, double elval, double ehval );

   /**
      add the bin width data, a pointer to an array with the bin upper edge information.
      This is needed when fitting with integral options
      The information is added for the previously inserted point.
      BinData::Add  must be called before
   */
   void AddBinUpEdge( const double* xup );

   /**
      return the value for the given fit point
   */
   double Value( unsigned int ipoint ) const
   {
      assert( ipoint < fMaxPoints );
      assert( fDataPtr );
      assert( fData.empty() || &fData.front() == fDataPtr );

      return fDataPtr[ipoint];
   }

   /**
      return a pointer to the value for the given fit point
   */
   const double *ValuePtr( unsigned int ipoint ) const
   {
      return &fDataPtr[ipoint];
   }

   /**
      return error on the value for the given fit point
      Safe (but slower) method returning correctly the error on the value
      in case of asymm errors return the average 0.5(eu + el)
   */

   const double * ErrorPtr(unsigned int ipoint) const{
      assert( ipoint < fMaxPoints );
      assert( kValueError == fErrorType || kCoordError == fErrorType ||
              kAsymError == fErrorType || kNoError == fErrorType );

      if ( fErrorType == kNoError )
         return nullptr;
      // assert( fErrorType == kCoordError );
      return &fDataErrorPtr[ ipoint ];
   }

   double Error( unsigned int ipoint ) const
   {
      assert( ipoint < fMaxPoints );
      assert( kValueError == fErrorType || kCoordError == fErrorType ||
              kAsymError == fErrorType || kNoError == fErrorType );

      if ( fErrorType == kNoError )
      {
         assert( !fDataErrorPtr && !fDataErrorHighPtr && !fDataErrorLowPtr );
         assert( fDataError.empty() && fDataErrorHigh.empty() && fDataErrorLow.empty() );
         return 1.0;
      }

      if ( fErrorType == kValueError ) // need to invert (inverror is stored)
      {
         assert( fDataErrorPtr && !fDataErrorHighPtr && !fDataErrorLowPtr );
         assert( fDataErrorHigh.empty() && fDataErrorLow.empty() );
         assert( fDataError.empty() || &fDataError.front() == fDataErrorPtr );

         double eval = fDataErrorPtr[ ipoint ];

         if (fWrapped)
            return eval;
         else
            return (eval != 0.0) ? 1.0/eval : 0.0;
      }

      if ( fErrorType == kAsymError )
      {  // return 1/2(el + eh)
         assert( !fDataErrorPtr && fDataErrorHighPtr && fDataErrorLowPtr );
         assert( fDataError.empty() );
         assert( fDataErrorHigh.empty() || &fDataErrorHigh.front() == fDataErrorHighPtr );
         assert( fDataErrorLow.empty() || &fDataErrorLow.front() == fDataErrorLowPtr );
         assert( fDataErrorLow.empty() == fDataErrorHigh.empty() );

         double eh = fDataErrorHighPtr[ ipoint ];
         double el = fDataErrorLowPtr[ ipoint ];

         return (el+eh) / 2.0;
      }

      assert( fErrorType == kCoordError );
      return fDataErrorPtr[ ipoint ];
   }

   void GetAsymError( unsigned int ipoint, double& lowError, double& highError ) const
   {
      assert( fErrorType == kAsymError );
      assert( !fDataErrorPtr && fDataErrorHighPtr && fDataErrorLowPtr );
      assert( fDataError.empty() );
      assert( fDataErrorHigh.empty() || &fDataErrorHigh.front() == fDataErrorHighPtr );
      assert( fDataErrorLow.empty() || &fDataErrorLow.front() == fDataErrorLowPtr );
      assert( fDataErrorLow.empty() == fDataErrorHigh.empty() );

      lowError = fDataErrorLowPtr[ ipoint ];
      highError = fDataErrorHighPtr[ ipoint ];
   }

   /**
      Return the inverse of error on the value for the given fit point
      useful when error in the coordinates are not stored and then this is used directly this as the weight in
      the least square function
   */
   double InvError( unsigned int ipoint ) const
   {
      assert( ipoint < fMaxPoints );
      assert( kValueError == fErrorType || kCoordError == fErrorType ||
              kAsymError == fErrorType || kNoError == fErrorType );

      if ( fErrorType == kNoError )
      {
         assert( !fDataErrorPtr && !fDataErrorHighPtr && !fDataErrorLowPtr );
         assert( fDataError.empty() && fDataErrorHigh.empty() && fDataErrorLow.empty() );
         return 1.0;
      }

      if ( fErrorType == kValueError ) // need to invert (inverror is stored)
      {
         assert( fDataErrorPtr && !fDataErrorHighPtr && !fDataErrorLowPtr );
         assert( fDataErrorHigh.empty() && fDataErrorLow.empty() );
         assert( fDataError.empty() || &fDataError.front() == fDataErrorPtr );

         double eval = fDataErrorPtr[ ipoint ];

         // in case of wrapped data the pointer stores the error and
         // not the inverse
         if (fWrapped)
            return 1.0 / eval;
         else
            return (eval != 0.0) ? eval : 0.0;
      }

      if ( fErrorType == kAsymError ) {
         // return inverse of 1/2(el + eh)
         assert( !fDataErrorPtr && fDataErrorHighPtr && fDataErrorLowPtr );
         assert( fDataError.empty() );
         assert( fDataErrorHigh.empty() || &fDataErrorHigh.front() == fDataErrorHighPtr );
         assert( fDataErrorLow.empty() || &fDataErrorLow.front() == fDataErrorLowPtr );
         assert( fDataErrorLow.empty() == fDataErrorHigh.empty() );

         double eh = fDataErrorHighPtr[ ipoint ];
         double el = fDataErrorLowPtr[ ipoint ];

         return 2.0 / (el+eh);
      }

      assert( fErrorType == kCoordError );
      // for coordinate error we store the error and not the inverse
      return 1.0 / fDataErrorPtr[ ipoint ];
   }


   /**
      retrieve at the same time a  pointer to the coordinate data and the fit value
      More efficient than calling Coords(i) and Value(i)
   */
   // not threadsafe, to be replaced with never constructs!
   // for example: just return std::array or std::vector, there's
   // is going to be only minor overhead in c++11.
   const double * GetPoint( unsigned int ipoint, double & value ) const
   {
      assert( ipoint < fMaxPoints );
      value = Value( ipoint );

      return Coords( ipoint );
   }

   /**
      returns a single coordinate error component of a point.
      This function is threadsafe in contrast to Coords(...)
      and can easily get vectorized by the compiler in loops
      running over the ipoint-index.
   */
   double GetCoordErrorComponent( unsigned int ipoint, unsigned int icoord ) const
   {
      assert( ipoint < fMaxPoints );
      assert( icoord < fDim );
      assert( fCoordErrorsPtr.size() == fDim );
      assert( fCoordErrorsPtr[icoord] );
      assert( fCoordErrors.empty() || &fCoordErrors[icoord].front() == fCoordErrorsPtr[icoord] );

      return fCoordErrorsPtr[icoord][ipoint];
   }

   /**
      Return a pointer to the errors in the coordinates for the given fit point
   */
   // not threadsafe, to be replaced with never constructs!
   // for example: just return std::array or std::vector, there's
   // is going to be only minor overhead in c++11.
   const double* CoordErrors( unsigned int ipoint ) const
   {
      assert( ipoint < fMaxPoints );
      assert( fpTmpCoordErrorVector );
      assert( fErrorType == kCoordError || fErrorType == kAsymError );

      for ( unsigned int i=0; i < fDim; i++ )
      {
         assert( fCoordErrorsPtr[i] );
         assert( fCoordErrors.empty() || &fCoordErrors[i].front() == fCoordErrorsPtr[i] );

         fpTmpCoordErrorVector[i] = fCoordErrorsPtr[i][ipoint];
      }

      return fpTmpCoordErrorVector;
   }


   /**
      retrieve in a single call a pointer to the coordinate data, value and inverse error for
      the given fit point.
      To be used only when type is kValueError or kNoError. In the last case the value 1 is returned
      for the error.
   */
   // not threadsafe, to be replaced with never constructs!
   // for example: just return std::array or std::vector, there's
   // is going to be only minor overhead in c++11.
   const double* GetPoint( unsigned int ipoint, double & value, double & invError ) const
   {
      assert( ipoint < fMaxPoints );
      assert( fErrorType == kNoError || fErrorType == kValueError );

      double e = Error( ipoint );

      if (fWrapped)
         invError = e;
      else
         invError = ( e != 0.0 ) ? 1.0/e : 1.0;

      return GetPoint( ipoint, value );
   }

   /**
      Retrieve the errors on the point (coordinate and value) for the given fit point
      It must be called only when the coordinate errors are stored otherwise it will produce an
      assert.
   */
   // not threadsafe, to be replaced with never constructs!
   // for example: just return std::array or std::vector, there's
   // is going to be only minor overhead in c++11.
   const double* GetPointError(unsigned int ipoint, double & errvalue) const
   {
      assert( ipoint < fMaxPoints );
      assert( fErrorType == kCoordError || fErrorType == kAsymError );

      errvalue = Error( ipoint );
      return CoordErrors( ipoint );
   }

   /**
      Get errors on the point (coordinate errors and asymmetric value errors) for the
      given fit point.
      It must be called only when the coordinate errors and asymmetric errors are stored
      otherwise it will produce an assert.
   */
   // not threadsafe, to be replaced with never constructs!
   // for example: just return std::array or std::vector, there's
   // is going to be only minor overhead in c++11.
   const double* GetPointError(unsigned int ipoint, double & errlow, double & errhigh) const
   {
      assert( ipoint < fMaxPoints );
      assert( fErrorType == kAsymError );
      assert( !fDataErrorPtr && fDataErrorHighPtr && fDataErrorLowPtr );
      assert( fDataError.empty() );
      assert( fDataErrorHigh.empty() || &fDataErrorHigh.front() == fDataErrorHighPtr );
      assert( fDataErrorLow.empty() || &fDataErrorLow.front() == fDataErrorLowPtr );
      assert( fDataErrorLow.empty() == fDataErrorHigh.empty() );

      errhigh = fDataErrorHighPtr[ ipoint ];
      errlow = fDataErrorLowPtr[ ipoint ];

      return CoordErrors( ipoint );
   }

   /**
      returns a single coordinate error component of a point.
      This function is threadsafe in contrast to Coords(...)
      and can easily get vectorized by the compiler in loops
      running over the ipoint-index.
   */
   double GetBinUpEdgeComponent( unsigned int ipoint, unsigned int icoord ) const
   {
      assert( icoord < fDim );
      assert( !fBinEdge.empty() );
      assert( ipoint < fBinEdge.front().size() );

      return fBinEdge[icoord][ipoint];
   }

   /**
      return an array containing the upper edge of the bin for coordinate i
      In case of empty bin they could be merged in a single larger bin
      Return a NULL pointer  if the bin width  is not stored
   */
   // not threadsafe, to be replaced with never constructs!
   // for example: just return std::array or std::vector, there's
   // is going to be only minor overhead in c++11.
   const double* BinUpEdge( unsigned int ipoint ) const
   {
      if ( fBinEdge.empty() || ipoint > fBinEdge.front().size() )
         return nullptr;

      GetBinUpEdgeCoordinates(ipoint, fpTmpBinEdgeVector);

      return fpTmpBinEdgeVector;
   }

   /**
    * Thread save version of function retrieving the bin up-edge in case of multidimensions
   */
   void GetBinUpEdgeCoordinates(unsigned int ipoint, double * x) const
   {
      if (fBinEdge.empty() || ipoint > fBinEdge.front().size()) return;
      assert(!fBinEdge.empty());
      assert(ipoint < fMaxPoints);
      for (unsigned int i = 0; i < fDim; i++) {
         x[i] = fBinEdge[i][ipoint];
      }
   }

   /**
      query if the data store the bin edges instead of the center
   */
   bool HasBinEdges() const {
      return fBinEdge.size() == fDim && fBinEdge[0].size() > 0;
   }

   /**
      retrieve the reference volume used to normalize the data when the option bin volume is set
   */
   double RefVolume() const { return fRefVolume; }

   /**
      set the reference volume used to normalize the data when the option bin volume is set
   */
   void SetRefVolume(double value) { fRefVolume = value; }

   /**
      retrieve the errortype
   */
   ErrorType GetErrorType( ) const
   {
      return fErrorType;
   }

   /**
      compute the total sum of the data content
      (sum of weights in case of weighted data set)
   */
   double SumOfContent() const { return fSumContent; }

   /**
      compute the total sum of the error square
      (sum of weight square in case of a weighted data set)
   */
   double SumOfError2() const { return fSumError2;}

   /**
      return true if the data set is weighted
      We cannot compute ourselfs because sometimes errors are filled with 1
      instead of zero (as in ROOT::Fit::FillData )
    */
   bool IsWeighted() const {
      return fIsWeighted;
   }

protected:
   void InitDataVector ();

   void InitializeErrors();

   void InitBinEdge();

   void UnWrap( );

   // compute sum of content and error squares
   void ComputeSums();

private:

   ErrorType fErrorType;
   bool fIsWeighted = false; // flag to indicate weighted data
   double fRefVolume;  // reference bin volume - used to normalize the bins in case of variable bins data
   double fSumContent = 0;  // total sum of the bin data content
   double fSumError2 = 0;  // total sum square of the errors

   /**
    * Stores the data values the same way as the coordinates.
    *
    */
   std::vector< double > fData;
   const double* fDataPtr;

   std::vector< std::vector< double > > fCoordErrors;
   std::vector< const double* > fCoordErrorsPtr;
   // This vector contains the coordinate errors
   // in the same way as fCoords.

   std::vector< double > fDataError;
   std::vector< double > fDataErrorHigh;
   std::vector< double > fDataErrorLow;
   const double*  fDataErrorPtr;
   const double*  fDataErrorHighPtr;
   const double*  fDataErrorLowPtr;
   // This vector contains the data error.
   // Either only fDataError or fDataErrorHigh and fDataErrorLow are used.

   double* fpTmpCoordErrorVector; // not threadsafe stuff!

   std::vector< std::vector< double > > fBinEdge;
   // vector containing the bin upper edge (coordinate will contain low edge)

   double* fpTmpBinEdgeVector; // not threadsafe stuff!
};


   } // end namespace Fit

} // end namespace ROOT



#endif /* ROOT_Fit_BinData */
