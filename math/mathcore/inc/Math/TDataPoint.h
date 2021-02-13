// @(#)root/mathcore:$Id$
// Authors: C. Gumpert    09/2011
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2011 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
//
// Header file for TDataPointclass
//

#ifndef ROOT_Math_TDataPoint
#define ROOT_Math_TDataPoint

//ROOT include(s)
#include "RtypesCore.h"


namespace ROOT
{
namespace Math
{


template<unsigned int K,typename _val_type = float>
class TDataPoint
{
public:
   typedef _val_type value_type;
   enum {
      kDimension = K //the dimensionality of this data point
   };
   static UInt_t Dimension() {return kDimension;}
   TDataPoint();
#ifndef __MAKECINT__
   template<typename _coord_typ>
   TDataPoint(const _coord_typ* pData,_val_type fWeight = 1);
#endif
   //virtual ~TDataPoint() {}
#ifndef __MAKECINT__
   template<typename _val>
   value_type   Distance(const TDataPoint<K,_val>& rPoint) const;
#endif
   value_type   GetCoordinate(unsigned int iAxis) const;
   value_type   GetWeight() const {return m_fWeight;}
   Bool_t       Less(TDataPoint& rPoint,unsigned int iAxis) const;
   void         SetCoordinate(unsigned int iAxis,_val_type fValue);
   void         SetWeight(float fWeight) {m_fWeight = fWeight;}

private:
   value_type   m_vCoordinates[K]; //coordinates
   value_type   m_fWeight;          //weight at this point
};

// some typedef definitions
typedef TDataPoint<1,Float_t>  TDataPoint1F;
typedef TDataPoint<2,Float_t>  TDataPoint2F;
typedef TDataPoint<3,Float_t>  TDataPoint3F;
typedef TDataPoint<1,Double_t> TDataPoint1D;
typedef TDataPoint<2,Double_t> TDataPoint2D;
typedef TDataPoint<3,Double_t> TDataPoint3D;

}//namespace Math
}//namespace ROOT

#include "Math/TDataPoint.icc"


#endif // ROOT_Math_TDataPoint
