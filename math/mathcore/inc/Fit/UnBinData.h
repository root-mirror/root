// @(#)root/mathcore:$Id$
// Author: L. Moneta Wed Aug 30 11:15:23 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class UnBinData

#ifndef ROOT_Fit_UnBinData
#define ROOT_Fit_UnBinData

#include "Fit/FitData.h"
#include "Math/Error.h"

#include <vector>

namespace ROOT {

   namespace Fit {


//___________________________________________________________________________________
/**
   Class describing the unbinned data sets (just x coordinates values) of any dimensions

              There is the option to construct UnBindata copying the data in (using the DataVector class)
              or using pointer to external data (DataWrapper) class.
              In general is found to be more efficient to copy the data.
              In case of really large data sets for limiting memory consumption then the other option can be used
              Specialized constructor exists for using external data up to 3 dimensions.

              When the data are copying in the number of points can be set later (or re-set) using Initialize and
              the data are inserted one by one using the Add method.
              It is mandatory to set the size before using the Add method.

             @ingroup  FitData
*/
class UnBinData : public FitData {

public :

  /**
    constructor from dimension of point  and max number of points (to pre-allocate vector)
  */

  explicit UnBinData( unsigned int maxpoints = 0, unsigned int dim = 1,
    bool isWeighted = false ) :
    FitData( maxpoints, isWeighted ? dim + 1 : dim ),
    fWeighted(isWeighted)
  {
    assert( dim >= 1 );
    assert( !fWeighted || dim >= 2 );
  }


  /**
    constructor from range and default option
  */
  explicit UnBinData ( const DataRange & range, unsigned int maxpoints = 0,
    unsigned int dim = 1, bool isWeighted = false ) :
    FitData( range, maxpoints, isWeighted ? dim + 1 : dim ),
    fWeighted(isWeighted)
  {
    assert( dim >= 1 );
    assert( !fWeighted || dim >= 2 );
  }

  /**
    constructor from options and range
  */
  UnBinData (const DataOptions & opt, const DataRange & range,
    unsigned int maxpoints = 0, unsigned int dim = 1, bool isWeighted = false ) :
    FitData( opt, range, maxpoints, isWeighted ? dim + 1 : dim ),
    fWeighted(isWeighted)
  {
    assert( dim >= 1 );
    assert( !fWeighted || dim >= 2 );
  }

  /**
    constructor for 1D external data (data are not copied inside)
  */
  UnBinData(unsigned int n, const double * dataX ) :
    FitData( n, dataX ),
    fWeighted( false )
  {
  }

  /**
    constructor for 2D external data (data are not copied inside)
    or 1D data with a weight (if isWeighted = true)
  */
  UnBinData(unsigned int n, const double * dataX, const double * dataY,
    bool isWeighted = false ) :
    FitData( n, dataX, dataY ),
    fWeighted( isWeighted )
  {
  }

  /**
    constructor for 3D external data (data are not copied inside)
    or 2D data with a weight (if isWeighted = true)
  */
  UnBinData(unsigned int n, const double * dataX, const double * dataY,
    const double * dataZ, bool isWeighted = false ) :
    FitData( n, dataX, dataY, dataZ ),
    fWeighted( isWeighted )
  {
  }

  /**
    constructor for multi-dim external data (data are not copied inside)
    Uses as argument an iterator of a list (or vector) containing the const double * of the data
    An example could be the std::vector<const double *>::begin
    In case of weighted data, the external data must have a dim+1 lists of data
    The apssed dim refers just to the coordinate size
  */
  template<class Iterator>
  UnBinData(unsigned int n, unsigned int dim, Iterator dataItr,
    bool isWeighted = false ) :
    FitData( n, isWeighted ? dim + 1 : dim, dataItr ),
    fWeighted( isWeighted )
  {
    assert( dim >= 1 );
    assert( !fWeighted || dim >= 2 );
  }

  /**
    constructor for 1D data and a range (data are copied inside according to the given range)
  */
  UnBinData(unsigned int maxpoints, const double * dataX, const DataRange & range) :
    FitData( range, maxpoints, dataX ),
    fWeighted( false )
  {
  }


  /**
    constructor for 2D data and a range (data are copied inside according to the given range)
    or 1 1D data set + weight. If is weighted  dataY is the pointer to the list of the weights
  */
  UnBinData(unsigned int maxpoints, const double * dataX, const double * dataY,
    const DataRange & range, bool isWeighted = false) :
    FitData( range, maxpoints, dataX, dataY ),
    fWeighted( isWeighted )
  {
  }

  /**
    constructor for 3D data and a range (data are copied inside according to the given range)
    or a 2D data set + weights. If is weighted  dataZ is the pointer to the list of the weights
  */
  UnBinData(unsigned int maxpoints, const double * dataX, const double * dataY,
    const double * dataZ, const DataRange & range, bool isWeighted = false) :
    FitData( range, maxpoints, dataX, dataY, dataZ ),
    fWeighted( isWeighted )
  {
  }

  /**
    constructor for multi-dim external data and a range (data are copied inside according to the range)
    Uses as argument an iterator of a list (or vector) containing the const double * of the data
    An example could be the std::vector<const double *>::begin
  */
  template<class Iterator>
  UnBinData( unsigned int maxpoints, unsigned int dim, Iterator dataItr, const DataRange & range, bool isWeighted = false ) :
    FitData( range, maxpoints, dim, dataItr ),
    fWeighted( isWeighted )
  {
  }

private:
  /// copy constructor (private)
  UnBinData(const UnBinData &) : FitData() { assert(false); }
  /// assignment operator  (private)
  UnBinData & operator= (const UnBinData &) { assert(false); return *this; }

public:
  /**
    destructor, delete pointer to internal data or external data wrapper
  */
  virtual ~UnBinData() {
  }

  /**
    preallocate a data set given size and dimension of the coordinates
    if a vector already exists with correct dimension (point size) extend the existing one
    to a total size of maxpoints (equivalent to a Resize)
  */
  //void Initialize(unsigned int maxpoints, unsigned int dim = 1, bool isWeighted = false);


  /**
    add one dim coordinate data (unweighted)
  */
  void Add(double x)
  {
    assert( !fWeighted );

    FitData::Add( x );
  }


  /**
    add 2-dim coordinate data
    can also be used to add 1-dim data with a weight
  */
  void Add(double x, double y)
  {
    assert( fDim == 2 );
    double dataTmp[] = { x, y };

    FitData::Add( dataTmp );
  }

  /**
    add 3-dim coordinate data
    can also be used to add 2-dim data with a weight
  */
  void Add(double x, double y, double z)
  {
    assert( fDim == 3 );
    double dataTmp[] = { x, y, z };

    FitData::Add( dataTmp );
  }

  /**
    add multi-dim coordinate data
  */
  void Add( const double* x )
  {
    FitData::Add( x );
  }

  /**
    add multi-dim coordinate data + weight
  */
  void Add(const double *x, double w)
  {
    assert( fWeighted );

    std::vector<double> tmpVec(fDim);
    std::copy( x, x + fDim - 1, tmpVec.begin() );
    tmpVec[fDim-1] = w;

    FitData::Add( &tmpVec.front() );
  }

  /**
    return weight
  */
  double Weight( unsigned int ipoint ) const
  {
    assert( ipoint < fNPoints );

    if ( !fWeighted ) return 1.0;
    return *GetCoordComponent(ipoint, fDim-1);
  }

  const double * WeightsPtr( unsigned int ipoint ) const
  {
    assert( ipoint < fNPoints );

    if ( !fWeighted ){
       MATH_ERROR_MSG("UnBinData::WeightsPtr","The function is unweighted!");
       return nullptr;
    }
    return GetCoordComponent(ipoint, fDim-1);
  }


  /**
    return coordinate data dimension
  */
  unsigned int NDim() const
  { return fWeighted ? fDim -1 : fDim; }

  bool IsWeighted() const
  {
    return fWeighted;
  }

  void Append( unsigned int newPoints, unsigned int dim = 1, bool isWeighted = false )
  {
    assert( !fWrapped );

    fWeighted = isWeighted;

    FitData::Append( newPoints, dim );
  }

private:
  bool fWeighted;

};


   } // end namespace Fit

} // end namespace ROOT



#endif /* ROOT_Fit_UnBinData */
