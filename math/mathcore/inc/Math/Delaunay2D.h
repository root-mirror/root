// @(#)root/mathcore:$Id: Delaunay2D.h,v 1.00
// Author: Daniel Funke, Lorenzo Moneta

/*************************************************************************
 * Copyright (C) 2015 ROOT Math Team                                     *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Header file for class Delaunay2D

#ifndef ROOT_Math_Delaunay2D
#define ROOT_Math_Delaunay2D

//for testing purposes HAS_CGAL can be [un]defined here
//#define HAS_CGAL

//for testing purposes THREAD_SAFE can [un]defined here
//#define THREAD_SAFE


#include "RtypesCore.h"

#include <map>
#include <vector>
#include <set>
#include <functional>

#ifdef HAS_CGAL
	/* CGAL uses the name PTR as member name in its Handle class
	 * but its a macro defined in mmalloc.h of ROOT
	 * Safe it, disable it and then re-enable it later on*/
	#pragma push_macro("PTR")
	#undef PTR

   #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
   #include <CGAL/Delaunay_triangulation_2.h>
   #include <CGAL/Triangulation_vertex_base_with_info_2.h>
   #include <CGAL/Interpolation_traits_2.h>
   #include <CGAL/natural_neighbor_coordinates_2.h>
   #include <CGAL/interpolation_functions.h>

	#pragma pop_macro("PTR")
#endif

#ifdef THREAD_SAFE
	#include<atomic> //atomic operations for thread safety
#endif


namespace ROOT {



   namespace Math {

/**

   Class to generate a Delaunay triangulation of a 2D set of points. 
   Algorithm based on **Triangle**, a two-dimensional quality mesh generator and 
   Delaunay triangulator from Jonathan Richard Shewchuk. 

   See [http://www.cs.cmu.edu/~quake/triangle.html]

   \ingroup MathCore
 */


class Delaunay2D  {

public:

	struct Triangle {
		double x[3];
		double y[3];
		UInt_t idx[3];

		#ifndef HAS_CGAL
		double invDenom; //see comment below in CGAL fall back section
		#endif
	};

	typedef std::vector<Triangle> Triangles;

public:

   
   Delaunay2D(int n, const double *x, const double * y, const double * z, double xmin=0, double xmax=0, double ymin=0, double ymax=0);

   /// set the input points for building the graph
   void SetInputPoints(int n, const double *x, const double * y, const double * z, double xmin=0, double xmax=0, double ymin=0, double ymax=0);

   /// Return the Interpolated z value corresponding to the (x,y) point
   double  Interpolate(double x, double y);

   /// Find all triangles 
   void      FindAllTriangles();

   /// return the number of triangles 
   Int_t     NumberOfTriangles() const {return fNdt;}

   double  XMin() const {return fXNmin;}
   double  XMax() const {return fXNmax;}
   double  YMin() const {return fYNmin;}
   double  YMax() const {return fYNmax;}

   /// set z value to be returned for points  outside the region 
   void      SetZOuterValue(double z=0.) { fZout = z; }

   /// return the user defined Z-outer value
   double ZOuterValue() const { return fZout; } 

   // iterators on the found triangles 
   Triangles::const_iterator begin() const { return fTriangles.begin(); }
   Triangles::const_iterator end()  const { return fTriangles.end(); }


private:

   // internal methods

   
   inline double Linear_transform(double x, double offset, double factor){
	   return (x+offset)*factor;
   }

   /// internal function to normalize the points
   void DoNormalizePoints();

   /// internal function to find the triangle
   /// use Triangle or CGAL if flag is set 
   void DoFindTriangles();

   /// internal method to compute the interpolation
   double  DoInterpolateNormalized(double x, double y);


   
private:
   // class is not copyable
   Delaunay2D(const Delaunay2D&); // Not implemented
   Delaunay2D& operator=(const Delaunay2D&); // Not implemented

protected:

   //typedef std::function<double(double)>                  Transformer;

   Int_t       fNdt;         //!Number of Delaunay triangles found
   Int_t       fNpoints;     //!Number of data points 

   const double   *fX;           //!Pointer to X array (managed externally)
   const double   *fY;           //!Pointer to Y array 
   const double   *fZ;           //!Pointer to Z array 

   double    fXNmin;       //!Minimum value of fXN
   double    fXNmax;       //!Maximum value of fXN
   double    fYNmin;       //!Minimum value of fYN
   double    fYNmax;       //!Maximum value of fYN

   //Transformer xTransformer; //!transform x values to mapped space
   //Transformer yTransformer; //!transform y values to mapped space

   double    fOffsetX;      //!Normalization offset X
   double    fOffsetY;      //!Normalization offset Y

   double    fScaleFactorX; //!Normalization factor X
   double    fScaleFactorY; //!Normalization factor Y

   double    fZout;        //!Height for points lying outside the convex hull

#ifdef THREAD_SAFE

   enum class Initialization : char {UNINITIALIZED, INITIALIZING, INITIALIZED};
   std::atomic<Initialization> fInit; //!Indicate initialization state

#else
   Bool_t      fInit;        //!True if FindAllTriangels() has been performed
#endif


   Triangles   fTriangles;   //!Triangles of Triangulation

#ifdef HAS_CGAL

   //Functor class for accessing the function values/gradients
   	template< class PointWithInfoMap, typename ValueType >
   	struct Data_access : public std::unary_function< typename PointWithInfoMap::key_type,
   				 std::pair<ValueType, bool> >
   	{

   	  Data_access(const PointWithInfoMap& points, const ValueType * values)
   			  : _points(points), _values(values){};

   	  std::pair< ValueType, bool>
   	  operator()(const typename PointWithInfoMap::key_type& p) const {
   		typename PointWithInfoMap::const_iterator mit = _points.find(p);
   		if(mit!= _points.end())
   		  return std::make_pair(_values[mit->second], true);
   		return std::make_pair(ValueType(), false);
   	  };

   	  const PointWithInfoMap& _points;
   	  const ValueType * _values;
   	};

   	typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
   	typedef CGAL::Triangulation_vertex_base_with_info_2<uint, K> Vb;
   	typedef CGAL::Triangulation_data_structure_2<Vb>             Tds;
   	typedef CGAL::Delaunay_triangulation_2<K, Tds>               Delaunay;
   	typedef CGAL::Interpolation_traits_2<K>                      Traits;
   	typedef K::FT                                                Coord_type;
   	typedef K::Point_2                                           Point;
   	typedef std::map<Point, Vb::Info, K::Less_xy_2>              PointWithInfoMap;
   	typedef Data_access< PointWithInfoMap, double >            Value_access;

   Delaunay fCGALdelaunay; //! CGAL delaunay triangulation object
   PointWithInfoMap fNormalizedPoints; //! Normalized function values

#else // HAS_CGAL
   //fallback to triangle library

   /* Using barycentric coordinates for inTriangle test and interpolation
    *
    * Given triangle ABC and point P, P can be expressed by
    *
    * P.x = la * A.x + lb * B.x + lc * C.x
    * P.y = la * A.y + lb * B.y + lc * C.y
    *
    * with lc = 1 - la - lb
    *
    * P.x = la * A.x + lb * B.x + (1-la-lb) * C.x
    * P.y = la * A.y + lb * B.y + (1-la-lb) * C.y
    *
    * Rearranging yields
    *
    * la * (A.x - C.x) + lb * (B.x - C.x) = P.x - C.x
    * la * (A.y - C.y) + lb * (B.y - C.y) = P.y - C.y
    *
    * Thus
    *
    * la = ( (B.y - C.y)*(P.x - C.x) + (C.x - B.x)*(P.y - C.y) ) / ( (B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y) )
    * lb = ( (C.y - A.y)*(P.x - C.x) + (A.x - C.x)*(P.y - C.y) ) / ( (B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y) )
    * lc = 1 - la - lb
    *
    * We save the inverse denominator to speedup computation
    *
    * invDenom = 1 / ( (B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y) )
    *
    * P is in triangle (including edges if
    *
    * 0 <= [la, lb, lc] <= 1
    *
    * The interpolation of P.z is
    *
    * P.z = la * A.z + lb * B.z + lc * C.z
    *
    */

   std::vector<double> fXN; //! normalized X
   std::vector<double> fYN; //! normalized Y

   /* To speed up localisation of points a grid is layed over normalized space
    *
    * A reference to triangle ABC is added to _all_ grid cells that include ABC's bounding box
    */

   static const int fNCells = 25; //! number of cells to divide the normalized space
   double fXCellStep; //! inverse denominator to calculate X cell = fNCells / (fXNmax - fXNmin)
   double fYCellStep; //! inverse denominator to calculate X cell = fNCells / (fYNmax - fYNmin)
   std::set<UInt_t> fCells[(fNCells+1)*(fNCells+1)]; //! grid cells with containing triangles

   inline unsigned int Cell(UInt_t x, UInt_t y) const {
	   return x*(fNCells+1) + y;
   }

   inline int CellX(double x) const {
	   return (x - fXNmin) * fXCellStep;
   }

   inline int CellY(double y) const {
	   return (y - fYNmin) * fYCellStep;
   }

#endif //HAS_CGAL


};


} // namespace Math
} // namespace ROOT

      
#endif
