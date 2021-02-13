// @(#)root/mathcore:$Id$
// Author: L. Moneta Wed Aug 30 11:05:02 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class DataRange

#ifndef ROOT_Fit_DataRange
#define ROOT_Fit_DataRange

#include <vector>
#include <utility>

namespace ROOT {

   namespace Fit {


//___________________________________________________________________________________
/**
   class describing the range in the coordinates
   it supports  multiple range in a coordinate.
   The rnage dimension is the dimension of the coordinate, its size is
   the number of interval for each coordinate.
   Default range is -inf, inf
   Range can be modified with the add range method

   @ingroup FitData
*/
class DataRange {

public:

   typedef std::vector<std::pair<double,double> > RangeSet;
   typedef std::vector< RangeSet >   RangeIntervals;

   /**
      Default constructor (infinite range)
   */
   explicit DataRange (unsigned int dim = 1) :
      fRanges ( std::vector<RangeSet> (dim) )
   {}

   /**
      construct a range for [xmin, xmax]
    */
   DataRange(double xmin, double xmax);

   /**
      construct a range for [xmin, xmax] , [ymin, ymax]
    */
   DataRange(double xmin, double xmax, double ymin, double ymax);
   /**
      construct a range for [xmin, xmax] , [ymin, ymax] , [zmin, zmax]
    */
   DataRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
   /**
      get range dimension
    */
   unsigned int NDim() const { return fRanges.size(); }

   /**
      return range size for coordinate icoord (starts from zero)
      Size == 0 indicates no range is present [-inf, + inf]
   */
   unsigned int Size(unsigned int icoord = 0) const {
      return icoord <  fRanges.size() ? fRanges[icoord].size() : 0;
   }

   /**
      return true if a range has been set in any of  the coordinates
      i.e. when  it is not [-inf,+inf] for all coordinates
      Avoid in case of multi-dim to loop on all the coordinated and ask the size
    */
   bool IsSet() const {
      for (unsigned int icoord = 0; icoord < fRanges.size(); ++icoord)
         if (fRanges[icoord].size() > 0) return true;
      return false;
   }

   /**
       return the vector of ranges for the coordinate icoord
   */
   const RangeSet & Ranges(unsigned int icoord = 0) const {
      // return icoord <  fRanges.size() ? fRanges[icoord] : RangeSet();
      return fRanges.at(icoord);
   }

   /**
       return the i-th range for the coordinate icoord.
       Useful method when only one range is present for the given coordinate
   */
   std::pair<double, double> operator() (unsigned int icoord = 0,unsigned int irange = 0) const;

   /**
      get the i-th range for given coordinate. If range does not exist
      return -inf, +inf
    */
   void GetRange(unsigned int irange, unsigned int icoord, double & xmin, double & xmax) const {
      if (Size(icoord)<= irange) GetInfRange(xmin,xmax);
      else {
         xmin = fRanges[icoord][irange].first;
         xmax = fRanges[icoord][irange].second;
      }
   }
   /**
      get the first range for given coordinate. If range does not exist
      return -inf, +inf
    */
   void GetRange(unsigned int icoord, double & xmin, double & xmax) const {
      if (Size(icoord) == 0) GetInfRange(xmin,xmax);
      else {
         xmin = fRanges[icoord].front().first;
         xmax = fRanges[icoord].front().second;
      }
   }
   /**
      get first range for the x - coordinate
    */
   void GetRange(double & xmin, double & xmax,unsigned int irange = 0) const {  GetRange(irange,0,xmin,xmax); }
   /**
      get range for the x and y coordinates
    */
   void GetRange(double & xmin, double & xmax, double & ymin, double & ymax,unsigned int irange = 0) const {
      GetRange(irange,0,xmin,xmax); GetRange(irange,1,ymin,ymax);
   }
   /**
      get range for the x and y and z coordinates
    */
   void GetRange(double & xmin, double & xmax, double & ymin, double & ymax, double & zmin, double & zmax,unsigned int irange=0) const {
      GetRange(irange,0,xmin,xmax); GetRange(irange,1,ymin,ymax); GetRange(irange,2,zmin,zmax);
   }
   /**
      get range for coordinates and fill the vector
    */
   void GetRange(double * xmin, double * xmax, unsigned int irange = 0)   const {
      for (unsigned int i = 0; i < fRanges.size(); ++i)
         GetRange(irange,i,xmin[i],xmax[i]);
   }

   /**
      Destructor (no operations)
   */
   ~DataRange ()  {}



   /**
      add a range [xmin,xmax] for the new coordinate icoord
      Adding a range does not delete existing one, but takes the OR with
      existing ranges.
      if want to replace range use method SetRange, which replace range with existing one
    */
   void AddRange(unsigned  int  icoord , double xmin, double xmax );

   /**
      add a range [xmin,xmax] for the first coordinate icoord
    */
   void AddRange(double xmin, double xmax ) { AddRange(0,xmin,xmax); }
   /**
      add a range [xmin,xmax] for the first and [ymin,ymax] for the second coordinate
    */
   void AddRange(double xmin, double xmax, double ymin, double ymax ) { AddRange(0,xmin,xmax); AddRange(1,ymin,ymax); }
   /**
      add a range [xmin,xmax] for the first and [ymin,ymax] for the second coordinate and
      [zmin,zmax] for the third coordinate
    */
   void AddRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax ) {
      AddRange(0,xmin,xmax); AddRange(1,ymin,ymax); AddRange(2,zmin,zmax); }

   /**
      set a range [xmin,xmax] for the new coordinate icoord
      If more range exists for other coordinates, delete the existing one and use it the new one
      Use Add range if want to keep the union of the existing ranges
    */
   void SetRange(unsigned  int  icoord , double xmin, double xmax );

   /**
      set a range [xmin,xmax] for the first coordinate icoord
    */
   void SetRange(double xmin, double xmax ) { SetRange(0,xmin,xmax); }
   /**
      set a range [xmin,xmax] for the first and [ymin,ymax] for the second coordinate
    */
   void SetRange(double xmin, double xmax, double ymin, double ymax ) { SetRange(0,xmin,xmax); SetRange(1,ymin,ymax); }
   /**
      set a range [xmin,xmax] for the first and [ymin,ymax] for the second coordinate and
      [zmin,zmax] for the third coordinate
    */
   void SetRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax ) {
      SetRange(0,xmin,xmax); SetRange(1,ymin,ymax); SetRange(2,zmin,zmax); }

   /**
      clear all ranges in one coordinate (is now -inf, +inf)
    */
   void Clear (unsigned  int  icoord = 0 );

   /**
      check if a point is inside the range for the given coordinate
    */
   bool IsInside(double x, unsigned int icoord = 0) const;

   /**
      check if a multi-dimpoint is inside the range 
    */
   bool IsInside(const double *x) const {
      bool ret = true;
      for (unsigned int idim = 0; idim < fRanges.size(); ++idim) { 
         ret &= IsInside(x[idim],idim);
         if (!ret) return ret;
      }
      return ret; 
   }

protected:
   /**
       internal function to remove all the existing ranges between xmin and xmax
       called when a new range is inserted
   */
   void CleanRangeSet(unsigned int icoord, double xmin, double xmax);

   // get the full range (-inf, +inf)
   static void GetInfRange(double &x1, double &x2);

private:

   RangeIntervals fRanges;  // list of all ranges


};

   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_DataRange */
