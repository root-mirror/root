// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableRearrangeTransform                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      rearrangement of input variables                                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland           *
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

#ifndef ROOT_TMVA_VariableRearrangeTransform
#define ROOT_TMVA_VariableRearrangeTransform

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// VariableRearrangeTransform                                           //
//                                                                      //
// rearrangement of input variables                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/VariableTransformBase.h"

#include <vector>

namespace TMVA {

   class VariableRearrangeTransform : public VariableTransformBase {

   public:

      typedef std::vector<Float_t>       FloatVector;

      VariableRearrangeTransform( DataSetInfo& dsi );
      virtual ~VariableRearrangeTransform( void );

      void   Initialize();
      Bool_t PrepareTransformation (const std::vector<Event*>&);

      virtual const Event* Transform(const Event* const, Int_t cls ) const;
      virtual const Event* InverseTransform( const Event* const, Int_t cls ) const;

      void WriteTransformationToStream ( std::ostream& ) const {}
      void ReadTransformationFromStream( std::istream&, const TString& ) { SetCreated(); }

      virtual void AttachXMLTo(void* parent);
      virtual void ReadFromXML( void* trfnode );

      virtual void PrintTransformation( std::ostream & o );

      // writer of function code
      virtual void MakeFunction( std::ostream& fout, const TString& fncName, Int_t part, UInt_t trCounter, Int_t cls );

      // provides string vector giving explicit transformation
      std::vector<TString>* GetTransformationStrings( Int_t cls ) const;

   private:

      ClassDef(VariableRearrangeTransform,0); // Variable transformation: normalization
   };

} // namespace TMVA

#endif 
