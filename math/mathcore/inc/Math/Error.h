// @(#)root/mathcore:$Id$
// Authors: L. Moneta

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2007 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Math_Error
#define ROOT_Math_Error

#ifdef DEBUG
#ifndef WARNINGMSG
#define WARNINGMSG
#endif
#endif



/**
   Pre-processor macro to report messages
   which can be configured to use ROOT error or
   simply an std::iostream in case of stan-alone builds
*/

#ifdef MATHCORE_STANDALONE

// use std::iostream instead of ROOT

#include <iostream>
#include <string>

#ifndef ROOT_MATH_LOG
#define ROOT_MATH_OS std::cerr
#else
#define ROOT_MATH_LOG
#endif

// giving a location + string

#define MATH_INFO_MSG(loc,str)                   \
   ROOT_MATH_OS << "Info in ROOT::Math::" << loc << ">: "  << str \
       << std::endl;
#define MATH_WARN_MSG(loc,str)                      \
   ROOT_MATH_OS << "Warning in ROOT::Math::" << loc << ">: " << str \
       << std::endl;
#define MATH_ERROR_MSG(loc,str)                   \
   ROOT_MATH_OS << "Error in ROOT::Math::" << loc << ">: " << str \
       << std::endl;

// giving location +  a value

# define MATH_INFO_VAL(loc,x)                                           \
   ROOT_MATH_OS << "Info in <ROOT::Math::" << loc << ">: " << #x << " = " << (x) << std::endl;
# define MATH_WARN_VAL(loc,x)                                           \
   ROOT_MATH_OS << "Warning in ROOT::Math::" << loc << ">: " << #x << " = " << (x) << std::endl;
# define MATH_ERROR_VAL(loc,x)                                          \
   ROOT_MATH_OS << "Error in ROOT::Math::" << loc << ">: " << #x << " = " << (x) << std::endl;

// giving a location + string + value

# define MATH_INFO_MSGVAL(loc,str,x)                                    \
   ROOT_MATH_OS << "Info in <ROOT::Math::" << loc << ">: "  << str << "; " << #x << " = " << (x) << std::endl;
# define MATH_WARN_MSGVAL(loc,str,x)                                    \
   ROOT_MATH_OS << "Warning in ROOT::Math::" << loc << ">: " << str << "; " << #x << " = " << (x) << std::endl;
# define MATH_ERROR_MSGVAL(loc,str,x)                                   \
   ROOT_MATH_OS << "Error in ROOT::Math::" << loc << ">: " << str << "; " << #x << " = " << (x) << std::endl;


#else
// use ROOT error reporting system

#include "TError.h"
#include "Math/Util.h"

#define  MATH_INFO_MSG(loc,str)                 \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
      ::Info(sl.c_str(),"%s",str);}
#define  MATH_WARN_MSG(loc,str)                 \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
      ::Warning(sl.c_str(),"%s",str);}
#define  MATH_ERROR_MSG(loc,str)                \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
      ::Error(sl.c_str(),"%s",str);}

# define MATH_INFO_VAL(loc,x)                                           \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
    std::string str = std::string(#x) + std::string(" = ") + ::ROOT::Math::Util::ToString(x); \
   ::Info(sl.c_str(),"%s",str.c_str() );}
# define MATH_WARN_VAL(loc,x)                                           \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
    std::string str = std::string(#x) + std::string(" = ") + ::ROOT::Math::Util::ToString(x); \
   ::Warning(sl.c_str(),"%s",str.c_str() );}
# define MATH_ERROR_VAL(loc,x)                                          \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
    std::string str = std::string(#x) + std::string(" = ") + ::ROOT::Math::Util::ToString(x); \
   ::Error(sl.c_str(),"%s",str.c_str() );}


# define MATH_INFO_MSGVAL(loc,txt,x)                    \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
    std::string str = std::string(txt) + std::string("; ") + std::string(#x) + std::string(" = ") + ::ROOT::Math::Util::ToString(x); \
   ::Info(sl.c_str(),"%s",str.c_str() );}
# define MATH_WARN_MSGVAL(loc,txt,x)                      \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
    std::string str = std::string(txt) + std::string("; ") + std::string(#x) + std::string(" = ") + ::ROOT::Math::Util::ToString(x); \
   ::Warning(sl.c_str(),"%s",str.c_str() );}
# define MATH_ERROR_MSGVAL(loc,txt,x)                     \
   {std::string sl = "ROOT::Math::" + std::string(loc); \
    std::string str = std::string(txt) + std::string("; ") + std::string(#x) + std::string(" = ") + ::ROOT::Math::Util::ToString(x); \
   ::Error(sl.c_str(),"%s",str.c_str() );}



#endif


#endif  // ROOT_MATH_Error
