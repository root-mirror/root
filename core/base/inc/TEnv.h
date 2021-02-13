// @(#)root/base:$Id$
// Author: Fons Rademakers   22/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEnv
#define ROOT_TEnv


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEnv                                                                 //
//                                                                      //
// The TEnv class reads config files, by default named .rootrc. Three   //
// types of config files are read: global, user and local files. The    //
// global file is $ROOTSYS/etc/system<name> (or ROOTETCDIR/system<name>)//
// the user file is $HOME/<name> and the local file is ./<name>.        //
// By setting the shell variable ROOTENV_NO_HOME=1 the reading of       //
// the $HOME/<name> resource file will be skipped. This might be useful //
// in case the home directory resides on an automounted remote file     //
// system and one wants to avoid this file system from being mounted.   //
//                                                                      //
// The format of the .rootrc file is similar to the .Xdefaults format:  //
//                                                                      //
//   [+]<SystemName>.<RootName|ProgName>.<name>[(type)]:  <value>       //
//                                                                      //
// Where <SystemName> is either Unix, WinNT, MacOS or Vms,              //
// <RootName> the name as given in the TApplication ctor (or "RootApp"  //
// in case no explicit TApplication derived object was created),        //
// <ProgName> the current program name and <name> the resource name,    //
// with optionally a type specification. <value> can be either a        //
// string, an integer, a float/double or a boolean with the values      //
// TRUE, FALSE, ON, OFF, YES, NO, OK, NOT. Booleans will be returned as //
// an integer 0 or 1. The options [+] allows the concatenation of       //
// values to the same resouce name.                                     //
//                                                                      //
// E.g.:                                                                //
//                                                                      //
//   Unix.Rint.Root.DynamicPath: .:$ROOTSYS/lib:~/lib                   //
//   myapp.Root.Debug:  FALSE                                           //
//   TH.Root.Debug: YES                                                 //
//   *.Root.MemStat: 1                                                  //
//                                                                      //
// <SystemName> and <ProgName> or <RootName> may be the wildcard "*".   //
// A # in the first column starts comment line.                         //
//                                                                      //
// For the currently defined resources (and their default values) see   //
// $ROOTSYS/etc/system.rootrc.                                          //
//                                                                      //
// Note that the .rootrc config files contain the config for all ROOT   //
// based applications.                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TString.h"

class THashList;
class TEnv;
class TEnvParser;
class TReadEnvParser;
class TWriteEnvParser;

enum EEnvLevel {
   kEnvGlobal,
   kEnvUser,
   kEnvLocal,
   kEnvChange,
   kEnvAll
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEnvRec                                                              //
//                                                                      //
// Individual TEnv records.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TEnvRec : public TObject {

friend class  TEnv;
friend class  TEnvParser;
friend class  TReadEnvParser;
friend class  TWriteEnvParser;

private:
   TString     fName;       // env rec key name
   TString     fType;       // env rec type
   TString     fValue;      // env rec value
   EEnvLevel   fLevel;      // env rec level
   Bool_t      fModified;   // if env rec has been modified

   TEnvRec(const char *n, const char *v, const char *t, EEnvLevel l);
   Int_t    Compare(const TObject *obj) const;
   void     ChangeValue(const char *v, const char *t, EEnvLevel l,
                        Bool_t append = kFALSE, Bool_t ignoredup = kFALSE);
   TString  ExpandValue(const char *v);

public:
   TEnvRec(): fName(), fType(), fValue(), fLevel(kEnvAll), fModified(kTRUE) { }
   ~TEnvRec();
   const char *GetName() const { return fName; }
   const char *GetValue() const { return fValue; }
   const char *GetType() const { return fType; }
   EEnvLevel   GetLevel() const { return fLevel; }
   ULong_t     Hash() const { return fName.Hash(); }

   ClassDef(TEnvRec,2)  // Individual TEnv records
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEnv                                                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TEnv : public TObject {

private:
   THashList        *fTable;     // hash table containing env records
   TString           fRcName;    // resource file base name
   Bool_t            fIgnoreDup; // ignore duplicates, don't issue warning

   TEnv(const TEnv&) = delete;
   TEnv& operator=(const TEnv&) = delete;

   const char       *Getvalue(const char *name) const;

public:
   TEnv(const char *name="");
   virtual ~TEnv();

   THashList          *GetTable() const { return fTable; }
   Bool_t              Defined(const char *name) const
                                    { return Getvalue(name) != nullptr; }

   virtual const char *GetRcName() const { return fRcName; }
   virtual void        SetRcName(const char *name) { fRcName = name; }

   virtual Int_t       GetValue(const char *name, Int_t dflt) const;
   virtual Double_t    GetValue(const char *name, Double_t dflt) const;
   virtual const char *GetValue(const char *name, const char *dflt) const;

   virtual void        SetValue(const char *name, const char *value,
                                EEnvLevel level = kEnvChange,
                                const char *type = nullptr);
   virtual void        SetValue(const char *name, EEnvLevel level = kEnvChange);
   virtual void        SetValue(const char *name, Int_t value);
   virtual void        SetValue(const char *name, Double_t value);

   virtual TEnvRec    *Lookup(const char *n) const;
   virtual Int_t       ReadFile(const char *fname, EEnvLevel level);
   virtual Int_t       WriteFile(const char *fname, EEnvLevel level = kEnvAll);
   virtual void        Save();
   virtual void        SaveLevel(EEnvLevel level);
   virtual void        Print(Option_t *option="") const;
   virtual void        PrintEnv(EEnvLevel level = kEnvAll) const;
   Bool_t              IgnoreDuplicates(Bool_t ignore);

   ClassDef(TEnv,2)  // Handle ROOT configuration resources
};

R__EXTERN TEnv *gEnv;

#endif
