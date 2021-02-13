// @(#)root/cont:$Id$
// Author: Philippe Canal 20/08/2003

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef R__RSTL_H
#define R__RSTL_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RStl                                                                 //
//                                                                      //
// Use to manage the code that needs to be generated for the STL        //
// by rootcint.  This class is reserved for rootcint and is a           //
// singleton.                                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <set>

#include "TClingUtils.h"

namespace clang {
   class CXXRecordDecl;
}

namespace cling {
   class Interpreter;
}

namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}

namespace ROOT {
namespace Internal {

   class RStl {
   private:
      typedef std::set<ROOT::TMetaUtils::AnnotatedRecordDecl,ROOT::TMetaUtils::AnnotatedRecordDecl::CompareByName> list_t;
      list_t fList;

   public:
      static RStl& Instance();
      ~RStl() {};

      void GenerateTClassFor(const char *requestedName, const clang::CXXRecordDecl *stlClass, const cling::Interpreter &interp, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt);
      void GenerateTClassFor(const clang::QualType &type, const cling::Interpreter &interp, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt);
      void Print();
      void WriteClassInit(std::ostream &strm,
                          const cling::Interpreter &interp,
                          const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt,
                          const ROOT::TMetaUtils::RConstructorTypes&,
                          bool &needCollectionProxy,
                          void (*emitStreamerInfo)(const char*) );
//       void WriteStreamer(FILE *file,const clang::CXXRecordDecl *stlcl);
//       void WriteStreamer(FILE *file);

   private:
      RStl() : fList() {};
      RStl(const RStl&);
      RStl& operator=(const RStl&);
   };

}
}
#endif // R__RSTL_H
