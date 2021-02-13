// @(#)root/gdml:$Id$
// Author: Andrei Gheata 05/12/2018

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGDMLMatrix
\ingroup Geometry_gdml
  This class is used in the process of reading and writing the GDML "matrix" tag.
It represents a matrix with arbitrary number of rows and columns, storing elements
in double precision.
*/

#include "TGDMLMatrix.h"

#include <cassert>

ClassImp(TGDMLMatrix)

//_____________________________________________________________________________
TGDMLMatrix::TGDMLMatrix(const char *name, size_t rows,size_t cols)
   : TNamed(name, "")
{
// Constructor
   if ((rows <= 0) || (cols <= 0))
   {
     Fatal("TGDMLMatrix::TGDMLMatrix(rows,cols)", "Wrong number of rows/cols");
   }
   fNrows = rows;
   fNcols = cols;
   fNelem = rows * cols;
   fMatrix = new Double_t[fNelem];
}

//_____________________________________________________________________________
TGDMLMatrix::TGDMLMatrix(const TGDMLMatrix& rhs)
  : TNamed(rhs), fNelem(rhs.fNelem), fNrows(rhs.fNrows), fNcols(rhs.fNcols), fMatrix(nullptr)
{
// Copy constructor
   if (rhs.fMatrix)
   {
     fMatrix = new Double_t[fNelem];
     memcpy(fMatrix, rhs.fMatrix, fNelem * sizeof(Double_t));
   }
}

//_____________________________________________________________________________
TGDMLMatrix& TGDMLMatrix::operator=(const TGDMLMatrix& rhs)
{
// Assignment
   if (this == &rhs)  { return *this; }
   TNamed::operator=(rhs);
   fNrows = rhs.fNrows;
   fNcols = rhs.fNcols;
   fNelem = fNrows * fNcols;
   if (rhs.fMatrix)
   {
      delete [] fMatrix;
      fMatrix = new Double_t[fNelem];
      memcpy(fMatrix, rhs.fMatrix, fNelem * sizeof(Double_t));
   }
   return *this;
}

//_____________________________________________________________________________
void TGDMLMatrix::Set(size_t r, size_t c, Double_t a)
{
   assert(r < fNrows && c < fNcols);
   fMatrix[fNcols*r+c] = a;
}

//_____________________________________________________________________________
Double_t TGDMLMatrix::Get(size_t r, size_t c) const
{
   assert(r < fNrows && c < fNcols);
   return fMatrix[fNcols*r+c];
}

//_____________________________________________________________________________
void TGDMLMatrix::Print(Option_t *) const
{
// Print info about this matrix
   printf("*** matrix: %-20s coldim = %zu  rows = %zu\n", GetName(), fNcols, fNrows);
   if (!fTitle.IsNull()) {
      printf("   %s\n", fTitle.Data());
      return;
   }
   for (size_t row = 0; row < fNrows; ++row) {
         printf("   ");
      for (size_t col = 0; col < fNcols; ++col) {
         printf("%8.3g", Get(row, col));
      }
      printf("\n");
   }
}
