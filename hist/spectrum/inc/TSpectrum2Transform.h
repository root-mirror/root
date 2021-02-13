// @(#)root/spectrum:$Id$
// Author: Miroslav Morhac   25/09/06

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TSpectrum2Transform
#define ROOT_TSpectrum2Transform

#include "TObject.h"

class TSpectrum2Transform : public TObject {
protected:
   Int_t      fSizeX;                      ///< x length of transformed data
   Int_t      fSizeY;                      ///< y length of transformed data
   Int_t      fTransformType;              ///< type of transformation (Haar, Walsh, Cosine, Sine, Fourier, Hartley, Fourier-Walsh, Fourier-Haar, Walsh-Haar, Cosine-Walsh, Cosine-Haar, Sine-Walsh, Sine-Haar)
   Int_t      fDegree;                     ///< degree of mixed transform, applies only for Fourier-Walsh, Fourier-Haar, Walsh-Haar, Cosine-Walsh, Cosine-Haar, Sine-Walsh, Sine-Haar transforms
   Int_t      fDirection;                  ///< forward or inverse transform
   Int_t      fXmin;                       ///< first channel x of filtered or enhanced region
   Int_t      fXmax;                       ///< last channel x of filtered or enhanced region
   Int_t      fYmin;                       ///< first channel y of filtered or enhanced region
   Int_t      fYmax;                       ///< last channel y of filtered or enhanced region
   Double_t   fFilterCoeff;                ///< value set in the filtered region
   Double_t   fEnhanceCoeff;               ///< multiplication coefficient applied in enhanced region;
public:
   enum {
       kTransformHaar =0,
       kTransformWalsh =1,
       kTransformCos =2,
       kTransformSin =3,
       kTransformFourier =4,
       kTransformHartley =5,
       kTransformFourierWalsh =6,
       kTransformFourierHaar =7,
       kTransformWalshHaar =8,
       kTransformCosWalsh =9,
       kTransformCosHaar =10,
       kTransformSinWalsh =11,
       kTransformSinHaar =12,
       kTransformForward =0,
       kTransformInverse =1
   };
   TSpectrum2Transform();
   TSpectrum2Transform(Int_t sizeX, Int_t sizeY);
   virtual ~TSpectrum2Transform();

protected:
   void                BitReverse(Double_t *working_space,Int_t num);
   void                BitReverseHaar(Double_t *working_space,Int_t shift,Int_t num,Int_t start);
   void                FourCos2(Double_t **working_matrix,Double_t *working_vector,Int_t numx,Int_t numy,Int_t direction,Int_t type);
   void                Fourier(Double_t *working_space,Int_t num,Int_t hartley,Int_t direction,Int_t zt_clear);
   void                General2(Double_t **working_matrix,Double_t *working_vector,Int_t numx,Int_t numy,Int_t direction,Int_t type,Int_t degree);
   Int_t               GeneralExe(Double_t *working_space,Int_t zt_clear,Int_t num,Int_t degree,Int_t type);
   Int_t               GeneralInv(Double_t *working_space,Int_t num,Int_t degree,Int_t type);
   void                Haar(Double_t *working_space,Int_t num,Int_t direction);
   void                HaarWalsh2(Double_t **working_matrix,Double_t *working_vector,Int_t numx,Int_t numy,Int_t direction,Int_t type);
   void                Walsh(Double_t *working_space,Int_t num);

public:
   void                Enhance(const Double_t **fSource, Double_t **fDest);
   void                FilterZonal(const Double_t **fSource, Double_t **fDest);
   void                SetDirection(Int_t direction);
   void                SetEnhanceCoeff(Double_t enhanceCoeff);
   void                SetFilterCoeff(Double_t filterCoeff);
   void                SetRegion(Int_t xmin, Int_t xmax, Int_t ymin, Int_t ymax);
   void                SetTransformType(Int_t transType, Int_t degree);
   void                Transform(const Double_t **fSource, Double_t **fDest);

   ClassDef(TSpectrum2Transform,1)  //Spectrum2 Transformer, it calculates classic orthogonal 2D transforms
};

#endif
