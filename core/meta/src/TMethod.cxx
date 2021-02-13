// @(#)root/meta:$Id$
// Author: Rene Brun   09/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TMethod
 Each ROOT class (see TClass) has a linked list of methods.
 This class describes one single method (member function).
 The method info is obtained via the CINT api. See class TCling.

 The method information is used a.o. by the THml class and by the
 TTree class.
*/

#include "strtok.h"
#include "strlcpy.h"
#include "snprintf.h"
#include "TClass.h"
#include "TList.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TMethodCall.h"
#include "TInterpreter.h"
#include "Strlen.h"
#include "TDataMember.h"


ClassImp(TMethod);

////////////////////////////////////////////////////////////////////////////////
/// Default TMethod ctor. TMethods are constructed in TClass.
/// Comment strings are pre-parsed to find out whether the method is
/// a context-menu item.

TMethod::TMethod(MethodInfo_t *info, TClass *cl) : TFunction(info)
{
   fClass        = cl;
   fGetterMethod = 0;
   fSetterMethod = 0;
   fMenuItem     = kMenuNoMenu;

   if (fInfo) {
      SetMenuItem(gCling->MethodInfo_Title(fInfo));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TMethod::TMethod(const TMethod& orig) : TFunction(orig)
{
   fClass        = orig.fClass;
   fMenuItem     = orig.fMenuItem;
   fGetter       = orig.fGetter;
   fGetterMethod = 0;
   fSetterMethod = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TMethod& TMethod::operator=(const TMethod& rhs)
{
   if (this != &rhs) {
      TFunction::operator=(rhs);
      fClass        = rhs.fClass;
      fMenuItem     = rhs.fMenuItem;
      fGetter       = rhs.fGetter;
      if (fGetterMethod)
         delete fGetterMethod;
      fGetterMethod = 0;
      if (fSetterMethod)
         delete fSetterMethod;
      fSetterMethod = 0;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TMethod::~TMethod()
{
   delete fGetterMethod;
   delete fSetterMethod;
}

////////////////////////////////////////////////////////////////////////////////
/// Clone method.

TObject *TMethod::Clone(const char *newname) const
{
   TNamed *newobj = new TMethod(*this);
   if (newname && strlen(newname)) newobj->SetName(newname);
   return newobj;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a comment string from the class declaration.

const char *TMethod::GetCommentString()
{
   return fInfo ? gCling->MethodInfo_Title(fInfo) : "";
}


////////////////////////////////////////////////////////////////////////////////
/// Using the CINT method arg information create a complete signature string.

void TMethod::CreateSignature()
{
   TFunction::CreateSignature();

   if (Property() & kIsConstMethod) fSignature += " const";
}

////////////////////////////////////////////////////////////////////////////////
/// Tries to guess DataMember from comment string
/// and Method's name <==(only if 1 Argument!).
/// If more then one argument=> returns pointer to the last argument.
/// It also sets MethodArgs' pointers to point to specified data members.
///
/// The form of comment string defining arguments is:
/// void XXX(Int_t x1, Float_t y2) //*ARGS={x1=>fX1,y2=>fY2}
/// where fX1, fY2 are data fields in the same class.
/// ("pointers" to data members)

TDataMember *TMethod::FindDataMember()
{
   Char_t *argstring = (char*)strstr(GetCommentString(),"*ARGS={");

   // the following statement has been commented (Rene). Not needed
   // it was making troubles in BuildRealData for classes with protected
   // default constructors.
   // if (!(GetClass()->GetListOfRealData())) GetClass()->BuildRealData();

   if (argstring) {

      // if we found any argument-specifying hints  - parse it

      if (!fMethodArgs) return 0;

      Int_t nchs = strlen(argstring);    // workspace...
      char *argstr = new char[nchs+1];   // workspace...
      char *ptr1 = 0;
      char *tok  = 0;
      char *ptr2 = 0;
      Int_t i;

      strlcpy(argstr,argstring,nchs+1);       //let's move it to "workspace"  copy
      char *rest;
      ptr2 = R__STRTOK_R(argstr, "{}", &rest); // extract the data!
      if (ptr2 == 0) {
         Fatal("FindDataMember","Internal error found '*ARGS=\"' but not \"{}\" in %s",GetCommentString());
         delete [] argstr;
         return 0;
      }
      ptr2 = R__STRTOK_R((char *)0, "{}", &rest);

      //extract argument tokens//
      char *tokens[20];
      Int_t cnt       = 0;
      Int_t token_cnt = 0;
      do {
         ptr1 = R__STRTOK_R((char *)(cnt++ ? 0 : ptr2), ",;", &rest); // extract tokens
                                                                   // separated by , or ;
         if (ptr1) {
            Int_t nch = strlen(ptr1);
            tok = new char[nch+1];
            strlcpy(tok,ptr1,nch+1);
            tokens[token_cnt] = tok;            //store this token.
            token_cnt++;
         }
      } while (ptr1);

      //now let's  parse all argument tokens...
      TClass     *cl = 0;
      TMethodArg *a  = 0;
      TMethodArg *ar = 0;
      TDataMember *member = 0;

      for (i=0; i<token_cnt;i++) {
         cnt = 0;
         ptr1 = R__STRTOK_R(tokens[i], "=>", &rest);         // LeftHandedSide=methodarg
         ptr2 = R__STRTOK_R((char *)0, "=>", &rest);         // RightHandedSide-points to datamember

         //find the MethodArg
         a      = 0;
         ar     = 0;
         member = 0;
         TIter nextarg(fMethodArgs);     // iterate through all arguments.
         while ((ar = (TMethodArg*)nextarg())) {
            if (!strcmp(ptr1,ar->GetName())) {
               a = ar;
               break;
            }
         }

         //now find the data member
         cl = GetClass()->GetBaseDataMember(ptr2);
         if (cl) {
            member = cl->GetDataMember(ptr2);
            if (a) a->fDataMember = member; //SET THE APROPRIATE FIELD !!!
                                     //We can do it - friend decl. in MethodArg
         }
         delete [] tokens[i];
      }
      delete [] argstr;
      return member; // nothing else to do! We return a pointer to the last
                     // found data member

   // if not found in comment string - try to guess it from name!
   } else {
      if (fMethodArgs)
         if (fMethodArgs->GetSize() != 1) return 0;

      TMethodArg *a = 0;
      if (fMethodArgs) a = (TMethodArg*)(fMethodArgs->First());

      char dataname[67]    = "";
      char basename[64]    = "";
      const char *funcname = GetName();
      if ( strncmp(funcname,"Get",3) == 0 || strncmp(funcname,"Set",3) == 0 )
         snprintf(basename,64,"%s",funcname+3);
      else if ( strncmp(funcname,"Is",2) == 0 )
         snprintf(basename,64,"%s",funcname+2);
      else if (strncmp(funcname, "Has", 3) == 0)
         snprintf(basename,64,"%s", funcname+3);
      else
         return 0;

      snprintf(dataname,67,"f%s",basename);

      TClass *cl = GetClass()->GetBaseDataMember(dataname);
      if (cl) {
         TDataMember *member   = cl->GetDataMember(dataname);
         if (a) a->fDataMember = member;
         return member;
      } else {
         snprintf(dataname,67,"fIs%s",basename);  //in case of IsEditable()
                                                        //and fIsEditable
         cl = GetClass()->GetBaseDataMember(dataname);
         if (cl) {
            TDataMember *member = cl->GetDataMember(dataname);
            if (a) a->fDataMember = member;
            return member;
         }
      }
   }

   //if nothing found - return null -pointer:
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return call environment for the getter method in case this is a
/// *TOGGLE method (for the context menu).

TMethodCall *TMethod::GetterMethod()
{
   if (!fGetterMethod && fMenuItem == kMenuToggle && fGetter != "" && fClass) {
      fGetterMethod = new TMethodCall(fClass, Getter(), "");
   }
   return fGetterMethod;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if this function object is pointing to a currently
/// loaded function.  If a function is unloaded after the TMethod
/// is created, the TMethod will be set to be invalid.

Bool_t TMethod::IsValid()
{
   // Register the transaction when checking the validity of the object.
   if (!fInfo && UpdateInterpreterStateMarker()) {
      DeclId_t newId = gInterpreter->GetFunction(fClass->GetClassInfo(), fName);
      if (newId) {
         MethodInfo_t *info = gInterpreter->MethodInfo_Factory(newId);
         Update(info);
      }
      return newId != 0;
   }
   return fInfo != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return call environment for this method in case this is a
/// *TOGGLE method which takes a single boolean or integer argument.

TMethodCall *TMethod::SetterMethod()
{
   if (!fSetterMethod && fMenuItem == kMenuToggle && fClass) {
      fSetterMethod = new TMethodCall(fClass, GetName(), "1");
   }
   return fSetterMethod;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns methodarg list and additionally updates fDataMember in TMethod by
/// calling FindDataMember();

TList *TMethod::GetListOfMethodArgs()
{
   if (!fMethodArgs){
      TFunction::GetListOfMethodArgs();
      FindDataMember();
   }
   return fMethodArgs;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the menu item as prescribed in the doctstring.

void TMethod::SetMenuItem(const char *docstring)
{
   if (docstring && strstr(docstring, "*TOGGLE")) {
      fMenuItem = kMenuToggle;
      const char *s;
      if ((s = strstr(docstring, "*GETTER="))) {
         fGetter = s+8;
         fGetter = fGetter.Strip(TString::kBoth);
      }
   } else
      if (docstring && strstr(docstring, "*MENU"))
         fMenuItem = kMenuDialog;
      else
         if (docstring && strstr(docstring, "*SUBMENU"))
            fMenuItem = kMenuSubMenu;
         else
            fMenuItem = kMenuNoMenu;
}

////////////////////////////////////////////////////////////////////////////////
/// Update the TMethod to reflect the new info.
///
/// This can be used to implement unloading (info == 0) and then reloading
/// (info being the 'new' decl address).

Bool_t TMethod::Update(MethodInfo_t *info)
{
   if (TFunction::Update(info)) {
      delete fGetterMethod; fGetterMethod = 0;
      delete fSetterMethod; fSetterMethod = 0;
      if (fInfo) {
         SetMenuItem(gCling->MethodInfo_Title(fInfo));
      }
      return kTRUE;
   } else {
      return kFALSE;
   }
}
