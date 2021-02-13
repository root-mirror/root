// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableIdentityTransform                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Identity transform                                                        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
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

#ifndef ROOT_TMVA_VariableIdentityTransform
#define ROOT_TMVA_VariableIdentityTransform

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// VariableIdentityTransform                                            //
//                                                                      //
// Linear interpolation class                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/VariableTransformBase.h"

#include <vector>

namespace TMVA {

   class VariableIdentityTransform : public VariableTransformBase {

   public:

      VariableIdentityTransform( DataSetInfo& dsi );
      virtual ~VariableIdentityTransform( void ) {}

      void   Initialize();
      Bool_t PrepareTransformation (const std::vector<Event*>& );

      void WriteTransformationToStream ( std::ostream& ) const {}
      void ReadTransformationFromStream( std::istream&, const TString& ) { SetCreated(); }

      virtual void AttachXMLTo(void* parent);
      virtual void ReadFromXML( void* trfnode );

      virtual const Event* Transform(const Event* const, Int_t cls ) const;
      virtual const Event* InverseTransform(const Event* const ev, Int_t cls ) const { return Transform( ev, cls ); }

      // writer of function code
      virtual void MakeFunction(std::ostream& fout, const TString& fncName, Int_t part, UInt_t trCounter, Int_t cls );

      ClassDef(VariableIdentityTransform,0); // Variable transformation: identity
   };

} // namespace TMVA

#endif 
