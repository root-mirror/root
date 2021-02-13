// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveRenderData
#define ROOT7_REveRenderData

#include <ROOT/REveVector.hxx>

#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {

class REveRenderData
{
private:
   std::string        fRnrFunc;
   std::vector<float> fVertexBuffer;
   std::vector<float> fNormalBuffer;
   std::vector<int>   fIndexBuffer;
   std::vector<float> fMatrix;

public:
   // If Primitive_e is changed, change also definition in EveElements.js.

   enum Primitive_e { GL_POINTS = 0, GL_LINES, GL_LINE_LOOP, GL_LINE_STRIP, GL_TRIANGLES };

   REveRenderData() = default;
   REveRenderData(const std::string &func, int size_vert = 0, int size_norm = 0, int size_idx = 0);

   void Reserve(int size_vert = 0, int size_norm = 0, int size_idx = 0);

   void PushV(float x) { fVertexBuffer.emplace_back(x); }

   void PushV(float x, float y, float z)
   {
      PushV(x);
      PushV(y);
      PushV(z);
   }

   void PushV(const REveVectorF &v)
   {
      PushV(v.fX);
      PushV(v.fY);
      PushV(v.fZ);
   }

   void PushV(float *v, int len) { fVertexBuffer.insert(fVertexBuffer.end(), v, v + len); }

   void PushN(float x) { fNormalBuffer.emplace_back(x); }

   void PushN(float x, float y, float z)
   {
      PushN(x);
      PushN(y);
      PushN(z);
   }

   void PushN(const REveVectorF &v)
   {
      PushN(v.fX);
      PushN(v.fY);
      PushN(v.fZ);
   }

   void PushI(UInt_t i) { fIndexBuffer.emplace_back(i); }

   void PushI(UInt_t i, UInt_t j, UInt_t k)
   {
      PushI(i);
      PushI(j);
      PushI(k);
   }

   void PushI(UInt_t *v, int len) { fIndexBuffer.insert(fIndexBuffer.end(), v, v + len); }

   void PushI(std::vector<UInt_t> &v) { fIndexBuffer.insert(fIndexBuffer.end(), v.begin(), v.end()); }

   void SetMatrix(const double *arr);

   const std::string GetRnrFunc() const { return fRnrFunc; }

   int SizeV() const { return fVertexBuffer.size(); }
   int SizeN() const { return fNormalBuffer.size(); }
   int SizeI() const { return fIndexBuffer.size(); }
   int SizeT() const { return fMatrix.size(); }

   int GetBinarySize() { return (SizeV() + SizeN() + SizeT()) * sizeof(float) + SizeI() * sizeof(int); }

   int Write(char *msg, int maxlen);
};

} // namespace Experimental
} // namespace ROOT

#endif
