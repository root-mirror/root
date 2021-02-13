// @(#)root/foam:$Id$
// Author: S. Jadach <mailto:Stanislaw.jadach@ifj.edu.pl>, P.Sawicki <mailto:Pawel.Sawicki@ifj.edu.pl>

//_____________________________________________________________________________
//                                                                            //
// Auxiliary class TFoamVect of n-dimensional vector, with dynamic allocation //
// used for the cartesian geometry of the TFoam  cells                        //
//                                                                            //
//_____________________________________________________________________________

#include "TFoamVect.h"

#include <iostream>
#include <iomanip>

/** \class TFoamVect

Auxiliary class TFoamVect of n-dimensional vector, with dynamic allocation
used for the cartesian geometry of the TFoam cells

*/

ClassImp(TFoamVect);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for streamer

TFoamVect::TFoamVect()
{
   fDim    =0;
   fCoords =0;
}

////////////////////////////////////////////////////////////////////////////////
/// User constructor creating n-dimensional vector
/// and allocating dynamically array of components

TFoamVect::TFoamVect(Int_t n)
{
   Int_t i;
   fDim=n;
   fCoords = 0;
   if (n>0) {
      fCoords = new Double_t[fDim];
      if(gDebug) {
         if(fCoords == 0)
            Error("TFoamVect", "Constructor failed to allocate\n");
      }
      for (i=0; i<n; i++) *(fCoords+i)=0.0;
   }
   if(gDebug>=3) Info("TFoamVect", "USER CONSTRUCTOR TFoamVect(const Int_t)\n ");
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TFoamVect::TFoamVect(const TFoamVect &Vect): TObject(Vect)
{
   fDim = Vect.fDim; fCoords = 0;
   if(fDim > 0)  fCoords = new Double_t[fDim];

   if(gDebug) {
      if(fCoords == 0) {
         Error("TFoamVect", "Constructor failed to allocate fCoords\n");
      }
   }

   for(Int_t i=0; i<fDim; i++)
      fCoords[i] = Vect.fCoords[i];

}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TFoamVect::~TFoamVect()
{
   if(gDebug>=3) Info("TFoamVect"," DESTRUCTOR TFoamVect~ \n");
   delete [] fCoords; //  free(fCoords)
   fCoords=0;
}

////////////////////////////////////////////////////////////////////////////////
/// substitution operator

TFoamVect& TFoamVect::operator =(const TFoamVect& Vect)
{
   Int_t i;
   if (&Vect == this) return *this;
   if( Vect.fDim < 0 )
       Error("TFoamVect","operator= : invalid  dimensions : %d and %d \n ",fDim,Vect.fDim);
   if( fDim != Vect.fDim ) {  // cleanup
      delete [] fCoords;
      fCoords = new Double_t[Vect.fDim];
   }
   fDim=Vect.fDim;
   for(i=0; i<fDim; i++)
      fCoords[i] = Vect.fCoords[i];
   if(gDebug>=3)  Info("TFoamVect", "SUBSITUTE operator =\n ");
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// [] is for access to elements as in ordinary matrix like a[j]=b[j]
/// (Perhaps against some strict rules but rather practical.)
/// Range protection is built in, consequently for substitution
/// one should use rather use a=b than explicit loop!

Double_t &TFoamVect::operator[](Int_t n)
{
   if ((n<0) || (n>=fDim)) {
      Error( "TFoamVect","operator[], out of range \n");
   }
   return fCoords[n];
}

////////////////////////////////////////////////////////////////////////////////
/// unary multiplication operator *=

TFoamVect& TFoamVect::operator*=(const Double_t &x)
{
   for(Int_t i=0;i<fDim;i++)
      fCoords[i] = fCoords[i]*x;
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// unary addition operator +=; adding vector c*=x,

TFoamVect& TFoamVect::operator+=(const TFoamVect& Shift)
{
   if( fDim != Shift.fDim){
      Error( "TFoamVect","operator+, different dimensions= %d %d \n",fDim,Shift.fDim);
   }
   for(Int_t i=0;i<fDim;i++)
      fCoords[i] = fCoords[i]+Shift.fCoords[i];
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// unary subtraction operator -=

TFoamVect& TFoamVect::operator-=(const TFoamVect& Shift)
{
   if( fDim != Shift.fDim) {
      Error( "TFoamVect","operator+, different dimensions= %d %d \n",fDim,Shift.fDim);
   }
   for(Int_t i=0;i<fDim;i++)
      fCoords[i] = fCoords[i]-Shift.fCoords[i];
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// addition operator +; sum of 2 vectors: c=a+b, a=a+b,
/// NEVER USE IT, VERY SLOW!!!

TFoamVect TFoamVect::operator+(const TFoamVect &p2)
{
   TFoamVect temp(fDim);
   temp  = (*this);
   temp += p2;
   return temp;
}

////////////////////////////////////////////////////////////////////////////////
/// subtraction operator -; difference of 2 vectors; c=a-b, a=a-b,
/// NEVER USE IT, VERY SLOW!!!

TFoamVect TFoamVect::operator-(const TFoamVect &p2)
{
   TFoamVect temp(fDim);
   temp  = (*this);
   temp -= p2;
   return temp;
}

////////////////////////////////////////////////////////////////////////////////
/// Loading in ordinary double prec. vector, sometimes can be useful

TFoamVect& TFoamVect::operator =(Double_t Vect[])
{
   Int_t i;
   for(i=0; i<fDim; i++)
      fCoords[i] = Vect[i];
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Loading in double prec. number, sometimes can be useful

TFoamVect& TFoamVect::operator =(Double_t x)
{
   if(fCoords != 0) {
      for(Int_t i=0; i<fDim; i++)
         fCoords[i] = x;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Printout of all vector components on "std::cout"

void TFoamVect::Print(Option_t *option) const
{
   if(!option) Error("Print ", "No option set \n");
   Int_t i;
   Int_t pr = std::cout.precision(7);
   std::cout << "(";
   for(i=0; i<fDim-1; i++) std::cout  << std::setw(12) << *(fCoords+i) << ",";
   std::cout  << std::setw(12) << *(fCoords+fDim-1);
   std::cout << ")";
   std::cout.precision(pr);
}
