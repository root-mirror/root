/* @(#)root/base:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2014, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RtypesCore
#define ROOT_RtypesCore

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RtypesCore                                                           //
//                                                                      //
// Basic types used by ROOT and required by TInterpreter.               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <ROOT/RConfig.hxx>

#include "DllImport.h"

#ifndef R__LESS_INCLUDES
#include <cstddef> // size_t, NULL
#endif

//---- Tag used by rootcling to determine constructor used for I/O.

class TRootIOCtor;

//---- types -------------------------------------------------------------------

typedef char           Char_t;      //Signed Character 1 byte (char)
typedef unsigned char  UChar_t;     //Unsigned Character 1 byte (unsigned char)
typedef short          Short_t;     //Signed Short integer 2 bytes (short)
typedef unsigned short UShort_t;    //Unsigned Short integer 2 bytes (unsigned short)
#ifdef R__INT16
typedef long           Int_t;       //Signed integer 4 bytes
typedef unsigned long  UInt_t;      //Unsigned integer 4 bytes
#else
typedef int            Int_t;       //Signed integer 4 bytes (int)
typedef unsigned int   UInt_t;      //Unsigned integer 4 bytes (unsigned int)
#endif
#ifdef R__B64    // Note: Long_t and ULong_t are currently not portable types
typedef int            Seek_t;      //File pointer (int)
typedef long           Long_t;      //Signed long integer 8 bytes (long)
typedef unsigned long  ULong_t;     //Unsigned long integer 8 bytes (unsigned long)
#else
typedef int            Seek_t;      //File pointer (int)
typedef long           Long_t;      //Signed long integer 4 bytes (long)
typedef unsigned long  ULong_t;     //Unsigned long integer 4 bytes (unsigned long)
#endif
typedef float          Float_t;     //Float 4 bytes (float)
typedef float          Float16_t;   //Float 4 bytes written with a truncated mantissa
typedef double         Double_t;    //Double 8 bytes
typedef double         Double32_t;  //Double 8 bytes in memory, written as a 4 bytes float
typedef long double    LongDouble_t;//Long Double
typedef char           Text_t;      //General string (char)
typedef bool           Bool_t;      //Boolean (0=false, 1=true) (bool)
typedef unsigned char  Byte_t;      //Byte (8 bits) (unsigned char)
typedef short          Version_t;   //Class version identifier (short)
typedef const char     Option_t;    //Option string (const char)
typedef int            Ssiz_t;      //String size (int)
typedef float          Real_t;      //TVector and TMatrix element type (float)
#if defined(R__WIN32) && !defined(__CINT__)
typedef __int64          Long64_t;  //Portable signed long integer 8 bytes
typedef unsigned __int64 ULong64_t; //Portable unsigned long integer 8 bytes
#else
typedef long long          Long64_t; //Portable signed long integer 8 bytes
typedef unsigned long long ULong64_t;//Portable unsigned long integer 8 bytes
#endif
typedef double         Axis_t;      //Axis values type (double)
typedef double         Stat_t;      //Statistics type (double)

typedef short          Font_t;      //Font number (short)
typedef short          Style_t;     //Style number (short)
typedef short          Marker_t;    //Marker number (short)
typedef short          Width_t;     //Line width (short)
typedef short          Color_t;     //Color number (short)
typedef short          SCoord_t;    //Screen coordinates (short)
typedef double         Coord_t;     //Pad world coordinates (double)
typedef float          Angle_t;     //Graphics angle (float)
typedef float          Size_t;      //Attribute size (float)

//---- constants ---------------------------------------------------------------

const Bool_t    kTRUE        = true;
const Bool_t    kFALSE       = false;

const Int_t     kMaxUChar    = 256;
const Int_t     kMaxChar     = kMaxUChar >> 1;
const Int_t     kMinChar     = -kMaxChar - 1;

const Int_t     kMaxUShort   = 65534;
const Int_t     kMaxShort    = kMaxUShort >> 1;
const Int_t     kMinShort    = -kMaxShort - 1;

const UInt_t    kMaxUInt     = UInt_t(~0);
const Int_t     kMaxInt      = Int_t(kMaxUInt >> 1);
const Int_t     kMinInt      = -kMaxInt - 1;

const ULong_t   kMaxULong    = ULong_t(~0);
const Long_t    kMaxLong     = Long_t(kMaxULong >> 1);
const Long_t    kMinLong     = -kMaxLong - 1;

const ULong64_t kMaxULong64  = ULong64_t(~0LL);
const Long64_t  kMaxLong64   = Long64_t(kMaxULong64 >> 1);
const Long64_t  kMinLong64   = -kMaxLong64 - 1;

const ULong_t   kBitsPerByte = 8;
const Ssiz_t    kNPOS        = ~(Ssiz_t)0;

//---- debug global ------------------------------------------------------------

R__EXTERN Int_t gDebug;


#endif
