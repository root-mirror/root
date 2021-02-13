// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveQuadSet.h"

#include "TEveManager.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TVirtualViewer3D.h"

/** \class TEveQuadSet
\ingroup TEve
Supports various internal formats that result in rendering of a
set of planar (lines, rectangles, hexagons with shared normal) objects.

Names of internal structures and their variables use A, B and C as
names for coordinate value-holders. Typical assignment is A->X,
B->Y, C->Z but each render mode can override this convention and
impose y or x as a fixed (third or C) coordinate. Alphabetic order
is obeyed in this correspondence.

For quad modes the deltas are expected to be positive.
For line modes negative deltas are ok.
*/

ClassImp(TEveQuadSet);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveQuadSet::TEveQuadSet(const char* n, const char* t) :
   TEveDigitSet   (n, t),

   fQuadType  (kQT_Undef),
   fDefWidth  (1),
   fDefHeight (1),
   fDefCoord  (0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveQuadSet::TEveQuadSet(EQuadType_e quadType, Bool_t valIsCol, Int_t chunkSize,
                         const char* n, const char* t) :
   TEveDigitSet   (n, t),

   fQuadType  (kQT_Undef),
   fDefWidth  (1),
   fDefHeight (1),
   fDefCoord  (0)
{
   Reset(quadType, valIsCol, chunkSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Return size of given atom type.

Int_t TEveQuadSet::SizeofAtom(TEveQuadSet::EQuadType_e qt)
{
   static const TEveException eH("TEveQuadSet::SizeofAtom ");

   switch (qt) {
      case kQT_Undef:                return 0;
      case kQT_FreeQuad:             return sizeof(QFreeQuad_t);
      case kQT_RectangleXY:
      case kQT_RectangleXZ:
      case kQT_RectangleYZ:          return sizeof(QRect_t);
      case kQT_RectangleXYFixedDim:  return sizeof(QRectFixDim_t);
      case kQT_RectangleXYFixedZ:
      case kQT_RectangleXZFixedY:
      case kQT_RectangleYZFixedX:    return sizeof(QRectFixC_t);
      case kQT_RectangleXYFixedDimZ:
      case kQT_RectangleXZFixedDimY:
      case kQT_RectangleYZFixedDimX: return sizeof(QRectFixDimC_t);
      case kQT_LineXZFixedY:
      case kQT_LineXYFixedZ:         return sizeof(QLineFixC_t);
      case kQT_HexagonXY:
      case kQT_HexagonYX:            return sizeof(QHex_t);
      default:                      throw(eH + "unexpected atom type.");
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the quad-set and reset the basic parameters.

void TEveQuadSet::Reset(TEveQuadSet::EQuadType_e quadType, Bool_t valIsCol,
                        Int_t chunkSize)
{
   fQuadType     = quadType;
   fValueIsColor = valIsCol;
   fDefaultValue = valIsCol ? 0 : kMinInt;
   if (fOwnIds)
      ReleaseIds();
   fPlex.Reset(SizeofAtom(fQuadType), chunkSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a quad specified with 4 vertices.

void TEveQuadSet::AddQuad(Float_t verts[12])
{
   static const TEveException eH("TEveQuadSet::AddQuad ");

   if (fQuadType != kQT_FreeQuad)
      throw(eH + "expect free quad-type.");

   QFreeQuad_t* fq = (QFreeQuad_t*) NewDigit();
   if (verts != 0)
     memcpy(fq->fVertices, verts, sizeof(fq->fVertices));
}

////////////////////////////////////////////////////////////////////////////////
/// Add a quad with a and b coordinates. Defaults are applied for
/// c coordinate and sizes.

void TEveQuadSet::AddQuad(Float_t a, Float_t b)
{
   AddQuad(a, b, fDefCoord, fDefWidth, fDefHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a quad with a, b and c coordinates. Defaults are applied
/// for sizes.

void TEveQuadSet::AddQuad(Float_t a, Float_t b, Float_t c)
{
   AddQuad(a, b, c, fDefWidth, fDefHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a quad with a and b coordinates and sizes. Default is applied
/// for c coordinate.

void TEveQuadSet::AddQuad(Float_t a, Float_t b, Float_t w, Float_t h)
{
   AddQuad(a, b, fDefCoord, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a quad with a, b and c coordinates and sizes.

void TEveQuadSet::AddQuad(Float_t a, Float_t b, Float_t c, Float_t w, Float_t h)
{
   static const TEveException eH("TEveQuadSet::AddAAQuad ");

   QOrigin_t& fq = * (QOrigin_t*) NewDigit();
   fq.fA = a; fq.fB = b;
   switch (fQuadType)
   {
      case kQT_RectangleXY:
      case kQT_RectangleXZ:
      case kQT_RectangleYZ:
      {
         QRect_t& q = (QRect_t&) fq;
         q.fC = c; q.fW = w; q.fH = h;
         break;
      }

      case kQT_RectangleXYFixedDim:
      {
         QRectFixDim_t& q =  (QRectFixDim_t&) fq;
         q.fC = c;
         break;
      }

      case kQT_RectangleXYFixedZ:
      case kQT_RectangleXZFixedY:
      case kQT_RectangleYZFixedX:
      {
         QRectFixC_t& q = (QRectFixC_t&) fq;
         q.fW = w; q.fH = h;
         break;
      }

      case kQT_RectangleXYFixedDimZ:
      case kQT_RectangleXZFixedDimY:
      case kQT_RectangleYZFixedDimX:
      {
         break;
      }

      default:
         throw(eH + "expect axis-aligned quad-type.");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add a line with starting coordinates and displacements.

void TEveQuadSet::AddLine(Float_t a, Float_t b, Float_t w, Float_t h)
{
   static const TEveException eH("TEveQuadSet::AddLine ");

   QOrigin_t& fq = * (QOrigin_t*) NewDigit();
   fq.fA = a; fq.fB = b;
   switch (fQuadType)
   {
      case kQT_LineXZFixedY:
      case kQT_LineXYFixedZ: {
         QLineFixC_t& q = (QLineFixC_t&) fq;
         q.fDx = w; q.fDy = h;
         break;
      }
      default:
         throw(eH + "expect line quad-type.");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add a hexagon with given center (a,b,c) and radius.

void TEveQuadSet::AddHexagon(Float_t a, Float_t b, Float_t c, Float_t r)
{
   static const TEveException eH("TEveQuadSet::AddHexagon ");

   QOrigin_t& fq = * (QOrigin_t*) NewDigit();
   fq.fA = a; fq.fB = b;
   switch (fQuadType)
   {
      case kQT_HexagonXY:
      case kQT_HexagonYX: {
         QHex_t& q = (QHex_t&) fq;
         q.fC = c; q.fR = r;
         break;
      }
      default:
         throw eH + "expects hexagon quad-type.";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill bounding-box information. Virtual from TAttBBox.
/// If member 'TEveFrameBox* fFrame' is set, frame's corners are
/// used as bbox.

void TEveQuadSet::ComputeBBox()
{
   static const TEveException eH("TEveQuadSet::ComputeBBox ");

   if (fFrame != 0)
   {
      BBoxInit();
      Int_t    n    = fFrame->GetFrameSize() / 3;
      Float_t *bbps = fFrame->GetFramePoints();
      for (int i=0; i<n; ++i, bbps+=3)
         BBoxCheckPoint(bbps);
   }
   else
   {
      if(fPlex.Size() == 0) {
         BBoxZero();
         return;
      }

      BBoxInit();
      if (fQuadType == kQT_RectangleXYFixedZ    ||
          fQuadType == kQT_RectangleXYFixedDimZ)
      {
         fBBox[4] = fDefCoord;
         fBBox[5] = fDefCoord;
      }
      else if (fQuadType == kQT_RectangleXZFixedY    ||
               fQuadType == kQT_RectangleXZFixedDimY)
      {
         fBBox[2] = fDefCoord;
         fBBox[3] = fDefCoord;
      }
      else if (fQuadType == kQT_RectangleYZFixedX    ||
               fQuadType == kQT_RectangleYZFixedDimX)
      {
         fBBox[0] = fDefCoord;
         fBBox[1] = fDefCoord;
      }

      TEveChunkManager::iterator qi(fPlex);

      switch (fQuadType)
      {

         case kQT_FreeQuad:
         {
            while (qi.next()) {
               const Float_t* p =  ((QFreeQuad_t*) qi())->fVertices;
               BBoxCheckPoint(p); p += 3;
               BBoxCheckPoint(p); p += 3;
               BBoxCheckPoint(p); p += 3;
               BBoxCheckPoint(p);
            }
            break;
         }

         case kQT_RectangleXY:
         {
            while (qi.next()) {
               QRect_t& q = * (QRect_t*) qi();
               if(q.fA        < fBBox[0]) fBBox[0] = q.fA;
               if(q.fA + q.fW > fBBox[1]) fBBox[1] = q.fA + q.fW;
               if(q.fB        < fBBox[2]) fBBox[2] = q.fB;
               if(q.fB + q.fH > fBBox[3]) fBBox[3] = q.fB + q.fH;
               if(q.fC        < fBBox[4]) fBBox[4] = q.fC;
               if(q.fC        > fBBox[5]) fBBox[5] = q.fC;
            }
            break;
         }

         case kQT_RectangleXZ:
         {
            while (qi.next()) {
               QRect_t& q = * (QRect_t*) qi();
               if(q.fA        < fBBox[0]) fBBox[0] = q.fA;
               if(q.fA + q.fW > fBBox[1]) fBBox[1] = q.fA + q.fW;
               if(q.fB        < fBBox[4]) fBBox[4] = q.fB;
               if(q.fB + q.fH > fBBox[5]) fBBox[5] = q.fB + q.fH;
               if(q.fC        < fBBox[2]) fBBox[2] = q.fC;
               if(q.fC        > fBBox[3]) fBBox[3] = q.fC;
            }
            break;
         }

         case kQT_RectangleYZ:
         {
            while (qi.next()) {
               QRect_t& q = * (QRect_t*) qi();
               if(q.fA        < fBBox[2]) fBBox[2] = q.fA;
               if(q.fA + q.fW > fBBox[3]) fBBox[3] = q.fA + q.fW;
               if(q.fB        < fBBox[4]) fBBox[4] = q.fB;
               if(q.fB + q.fH > fBBox[5]) fBBox[5] = q.fB + q.fH;
               if(q.fC        < fBBox[0]) fBBox[0] = q.fC;
               if(q.fC        > fBBox[1]) fBBox[1] = q.fC;
            }
            break;
         }

         case kQT_RectangleXYFixedDim:
         {
            const Float_t& w = fDefWidth;
            const Float_t& h = fDefHeight;
            while (qi.next()) {
               QRectFixDim_t& q = * (QRectFixDim_t*) qi();
               if(q.fA     < fBBox[0]) fBBox[0] = q.fA;
               if(q.fA + w > fBBox[1]) fBBox[1] = q.fA + w;
               if(q.fB     < fBBox[2]) fBBox[2] = q.fB;
               if(q.fB + h > fBBox[3]) fBBox[3] = q.fB + h;
               if(q.fC     < fBBox[4]) fBBox[4] = q.fC;
               if(q.fC     > fBBox[5]) fBBox[5] = q.fC;
            }
            break;
         }

         case kQT_RectangleXYFixedZ:
         {
            while (qi.next()) {
               QRectFixC_t& q = * (QRectFixC_t*) qi();
               if(q.fA        < fBBox[0]) fBBox[0] = q.fA;
               if(q.fA + q.fW > fBBox[1]) fBBox[1] = q.fA + q.fW;
               if(q.fB        < fBBox[2]) fBBox[2] = q.fB;
               if(q.fB + q.fH > fBBox[3]) fBBox[3] = q.fB + q.fH;
            }
            break;
         }

         case kQT_RectangleXZFixedY:
         {
            while (qi.next()) {
               QRectFixC_t& q = * (QRectFixC_t*) qi();
               if(q.fA        < fBBox[0]) fBBox[0] = q.fA;
               if(q.fA + q.fW > fBBox[1]) fBBox[1] = q.fA + q.fW;
               if(q.fB        < fBBox[4]) fBBox[4] = q.fB;
               if(q.fB + q.fH > fBBox[5]) fBBox[5] = q.fB + q.fH;
            }
            break;
         }

         case kQT_RectangleYZFixedX:
         {
            while (qi.next()) {
               QRectFixC_t& q = * (QRectFixC_t*) qi();
               if(q.fA        < fBBox[2]) fBBox[2] = q.fA;
               if(q.fA + q.fW > fBBox[3]) fBBox[3] = q.fA + q.fW;
               if(q.fB        < fBBox[4]) fBBox[4] = q.fB;
               if(q.fB + q.fH > fBBox[5]) fBBox[5] = q.fB + q.fH;
            }
            break;
         }

         case kQT_RectangleXYFixedDimZ:
         {
            const Float_t& w = fDefWidth;
            const Float_t& h = fDefHeight;
            while (qi.next()) {
               QRectFixDimC_t& q = * (QRectFixDimC_t*) qi();
               if(q.fA     < fBBox[0]) fBBox[0] = q.fA;
               if(q.fA + w > fBBox[1]) fBBox[1] = q.fA + w;
               if(q.fB     < fBBox[2]) fBBox[2] = q.fB;
               if(q.fB + h > fBBox[3]) fBBox[3] = q.fB + h;
            }
            break;
         }

         case kQT_RectangleXZFixedDimY:
         {
            const Float_t& w = fDefWidth;
            const Float_t& h = fDefHeight;
            while (qi.next()) {
               QRectFixDimC_t& q = * (QRectFixDimC_t*) qi();
               if(q.fA     < fBBox[0]) fBBox[0] = q.fA;
               if(q.fA + w > fBBox[1]) fBBox[1] = q.fA + w;
               if(q.fB     < fBBox[4]) fBBox[4] = q.fB;
               if(q.fB + h > fBBox[5]) fBBox[5] = q.fB + h;
            }
            break;
         }

         case kQT_RectangleYZFixedDimX:
         {
            const Float_t& w = fDefWidth;
            const Float_t& h = fDefHeight;
            while (qi.next()) {
               QRectFixDimC_t& q = * (QRectFixDimC_t*) qi();
               if(q.fA     < fBBox[2]) fBBox[2] = q.fA;
               if(q.fA + w > fBBox[3]) fBBox[3] = q.fA + w;
               if(q.fB     < fBBox[4]) fBBox[4] = q.fB;
               if(q.fB + h > fBBox[5]) fBBox[5] = q.fB + h;
            }
            break;
         }

         // TEveLine modes

         case kQT_LineXYFixedZ:
         {
            while (qi.next()) {
               QLineFixC_t& q = * (QLineFixC_t*) qi();
               BBoxCheckPoint(q.fA,         q.fB,         fDefCoord);
               BBoxCheckPoint(q.fA + q.fDx, q.fB + q.fDy, fDefCoord);
            }
            break;
         }

         case kQT_LineXZFixedY:
         {
            while (qi.next()) {
               QLineFixC_t& q = * (QLineFixC_t*) qi();
               BBoxCheckPoint(q.fA,         fDefCoord, q.fB);
               BBoxCheckPoint(q.fA + q.fDx, fDefCoord, q.fB + q.fDy);
            }
            break;
         }

         // Hexagon modes

         // Ignore 'slight' difference, assume square box for both cases.
         case kQT_HexagonXY:
         case kQT_HexagonYX:
         {
            while (qi.next()) {
               QHex_t& q = * (QHex_t*) qi();
               BBoxCheckPoint(q.fA-q.fR, q.fB-q.fR, q.fC);
               BBoxCheckPoint(q.fA+q.fR, q.fB+q.fR, q.fC);
            }
            break;
         }

         default:
         {
            throw(eH + "unsupported quad-type.");
         }

      } // end switch quad-type
   } // end if frame ... else ...
}
