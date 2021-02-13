// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Volume                                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Volume for BinarySearchTree                                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_Volume
#define ROOT_TMVA_Volume

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Volume                                                               //
//                                                                      //
// Volume for BinarySearchTree                                          //
//                                                                      //
// volume element: cubic variable space beteen upper and lower bonds of //
// nvar-dimensional variable space                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include "RtypesCore.h"

namespace TMVA {
  
   class Volume {

   public:

      // constructors
      Volume( std::vector<Float_t>* l, std::vector<Float_t>* u = 0);
      Volume( std::vector<Double_t>* l = 0, std::vector<Double_t>* u = 0);
      Volume( Volume& );
      Volume( Float_t* l , Float_t* u , Int_t nvar );
      Volume( Double_t* l , Double_t* u , Int_t nvar );
      Volume( Float_t l , Float_t u );
      Volume( Double_t l , Double_t u );

      // destructor
      virtual ~Volume( void );

      // operators
      Volume& operator=( const Volume& );

      // destruct the volue 
      void Delete       ( void );
      // "scale" the volume by multiplying each upper and lower boundary by "f" 
      void Scale        ( Double_t f );
      // "scale" the volume by symmetrically blowing up the interval in each dimension
      void ScaleInterval( Double_t f );
      void Print        ( void ) const;

      // allow direct access for better speed
      std::vector<Double_t> *fLower;    // vector with lower volume dimensions
      std::vector<Double_t> *fUpper;    // vector with upper volume dimensions

   private:

      Bool_t                fOwnerShip; // flag if "boundary vector" is owned by the volume of not
   };

} // namespace TMVA

#endif
