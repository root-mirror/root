// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariablePCATransform                                                  *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Principal value composition of input variables                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland           *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
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

#ifndef ROOT_TMVA_VariablePCATransform
#define ROOT_TMVA_VariablePCATransform

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// VariablePCATransform                                                 //
//                                                                      //
// Linear interpolation class                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPrincipal.h"

#include "TMVA/VariableTransformBase.h"

#include <vector>

namespace TMVA {

   class VariablePCATransform : public VariableTransformBase {

   public:
  
      VariablePCATransform( DataSetInfo& dsi );
      virtual ~VariablePCATransform( void );

      void   Initialize();
      Bool_t PrepareTransformation (const std::vector<Event*>&);

      virtual const Event* Transform(const Event* const, Int_t cls ) const;
      virtual const Event* InverseTransform(const Event* const, Int_t cls ) const;

      void WriteTransformationToStream ( std::ostream& ) const;
      void ReadTransformationFromStream( std::istream&, const TString& );

      virtual void AttachXMLTo(void* parent);
      virtual void ReadFromXML( void* trfnode );

      // writer of function code
      virtual void MakeFunction( std::ostream& fout, const TString& fncName, Int_t part, UInt_t trCounter, Int_t cls );

   private:

      void CalculatePrincipalComponents( const std::vector< Event*>& );
      void X2P( std::vector<Float_t>&, const std::vector<Float_t>&, Int_t cls ) const;
      void P2X( std::vector<Float_t>&, const std::vector<Float_t>&, Int_t cls ) const;

      // store relevant parts of PCA locally
      std::vector<TVectorD*> fMeanValues;   // mean values
      std::vector<TMatrixD*> fEigenVectors; // eigenvectors

      ClassDef(VariablePCATransform,0); // Variable transformation: Principal Value Composition
   };

} // namespace TMVA

#endif 

