// @(#)root/tree:$Id$
// Author:  Lukasz Janyst <ljanyst@cern.ch> 23/01/2008

//------------------------------------------------------------------------------
// file:   TIndArray.h
//------------------------------------------------------------------------------

#ifndef ROOT_TIndArray
#define ROOT_TIndArray

#include "RtypesCore.h"


class TIndArray
{
   public:
      TIndArray():
         fElems( 0 ), fCapacity( 0 ), fArr( nullptr ) {};

      virtual ~TIndArray()
      {
         delete [] fArr;
      }

      void ClearAndResize( UInt_t size )
      {
         delete [] fArr;
         fElems = 0;
         fArr = new UChar_t[size];
         fCapacity = size;
      }

      UInt_t   GetCapacity() { return fCapacity; }
      UInt_t   GetNumItems() { return fElems; }
      void     SetNumItems( UInt_t items ) { fElems = items;}
      UChar_t &At( Int_t ind ) { return fArr[ind]; }
      void     Clear() { fElems = 0; }

   private:
      UInt_t   fElems;     ///<  Number of elements stored in the array
      UInt_t   fCapacity;  ///<! Capacity of the array
      UChar_t *fArr;       ///<[fElems] The array

};

#endif // ROOT_TIndArray
