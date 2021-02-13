// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableNormalizeTransform                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Decorrelation of input variables                                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
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

#ifndef ROOT_TMVA_VariableNormalizeTransform
#define ROOT_TMVA_VariableNormalizeTransform

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// VariableNormalizeTransform                                           //
//                                                                      //
// Linear interpolation class                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMatrixDfwd.h"

#include "TMVA/VariableTransformBase.h"

#include <vector>

namespace TMVA {

   class VariableNormalizeTransform : public VariableTransformBase {

   public:

      typedef std::vector<Float_t>       FloatVector;
      typedef std::vector< FloatVector > VectorOfFloatVectors;
      VariableNormalizeTransform( DataSetInfo& dsi, TString strcor=""  );
      virtual ~VariableNormalizeTransform( void );

      void   Initialize();
      Bool_t PrepareTransformation (const std::vector<Event*>&);

      virtual const Event* Transform(const Event* const, Int_t cls ) const;
      virtual const Event* InverseTransform( const Event* const, Int_t cls ) const;

      void WriteTransformationToStream ( std::ostream& ) const;
      void ReadTransformationFromStream( std::istream&, const TString& );
      void BuildTransformationFromVarInfo( const std::vector<TMVA::VariableInfo>& var );

      virtual void AttachXMLTo(void* parent);
      virtual void ReadFromXML( void* trfnode );

      virtual void PrintTransformation( std::ostream & o );

      // writer of function code
      virtual void MakeFunction( std::ostream& fout, const TString& fncName, Int_t part, UInt_t trCounter, Int_t cls );

      // provides string vector giving explicit transformation
      std::vector<TString>* GetTransformationStrings( Int_t cls ) const;

   private:

      Bool_t           fNoOffset;

      void CalcNormalizationParams( const std::vector< Event*>& events);

      //      mutable Event*           fTransformedEvent;

      VectorOfFloatVectors                   fMin;       //! Min of source range
      VectorOfFloatVectors                   fMax;       //! Max of source range

      ClassDef(VariableNormalizeTransform,0); // Variable transformation: normalization
   };

} // namespace TMVA

#endif 
