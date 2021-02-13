// @(#)root/tmva $Id$    
// Author: Andrzej Zemla

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : SVKernelMatrix                                                        *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Kernel matrix for Support Vector Machine                                  *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Marcin Wolter  <Marcin.Wolter@cern.ch> - IFJ PAN, Krakow, Poland          *
 *      Andrzej Zemla  <azemla@cern.ch>        - IFJ PAN, Krakow, Poland          *
 *      (IFJ PAN: Henryk Niewodniczanski Inst. Nucl. Physics, Krakow, Poland)     *   
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      PAN, Krakow, Poland                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_SVKernelMatrix
#define ROOT_TMVA_SVKernelMatrix

#include "RtypesCore.h"

#include <vector>

namespace TMVA {

   class SVEvent;
   class SVKernelFunction;
   class MsgLogger;

   class SVKernelMatrix {

   public:

      //constructors
      SVKernelMatrix();
      SVKernelMatrix( std::vector<TMVA::SVEvent*>*, SVKernelFunction* );
      
      //destructor
      ~SVKernelMatrix();
      
      //functions
      Float_t* GetLine   ( UInt_t );
      Float_t* GetColumn ( UInt_t col ) { return this->GetLine(col);}
      Float_t  GetElement( UInt_t i, UInt_t j );

   private:

      UInt_t               fSize;              // matrix size
      SVKernelFunction*    fKernelFunction;    // kernel function
      Float_t**            fSVKernelMatrix;    // kernel matrix

      mutable MsgLogger* fLogger;                     //! message logger
      MsgLogger& Log() const { return *fLogger; }

   };
}
#endif
