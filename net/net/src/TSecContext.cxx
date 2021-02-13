// @(#)root/net:$Id$
// Author: G. Ganis   19/03/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSecContext                                                          //
//                                                                      //
// Contains details about an established security context               //
// Used by THostAuth                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include <cstdlib>

#include "strlcpy.h"
#include "snprintf.h"
#include "TSecContext.h"
#include "TSocket.h"
#include "TUrl.h"
#include "TROOT.h"
#include "TError.h"
#include "TVirtualMutex.h"

ClassImp(TSecContext);
ClassImp(TSecContextCleanup);

const TDatime kROOTTZERO = 788914800;

////////////////////////////////////////////////////////////////////////////////
/// Ctor for SecContext object.

TSecContext::TSecContext(const char *user, const char *host, Int_t meth,
                         Int_t offset, const char *id,
                         const char *token, TDatime expdate, void *ctx)
            : TObject()
{
   R__ASSERT(gROOT);

   fContext = ctx;
   fCleanup = new TList;
   fExpDate = expdate;
   if (offset > -1) {
      if (fExpDate < TDatime()) {
         // This means expdate was not initialized
         // We set it to default, ie 1 day from now
         fExpDate.Set(TDatime().GetDate() + 1, TDatime().GetTime());
      }
   }
   fHost    = host;
   fID      = id;
   fMethod  = meth;
   fMethodName = "";
   fOffSet  = offset;
   fToken   = token;
   fUser    = user;

   // Keep official list updated with active TSecContexts
   if (fOffSet > -1) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSecContexts()->Add(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Ctor for SecContext object.
/// User and host from url = user@host .

TSecContext::TSecContext(const char *url, Int_t meth, Int_t offset,
                         const char *token, const char *id,
                         TDatime expdate, void *ctx)
            : TObject()
{
   R__ASSERT(gROOT);

   fContext = ctx;
   fCleanup = new TList;
   fExpDate = expdate;
   if (offset > -1) {
      if (fExpDate < TDatime()) {
         // This means expdate was not initialized
         // We set it to default, ie 1 day from now
         fExpDate.Set(TDatime().GetDate() + 1, TDatime().GetTime());
      }
   }
   fHost    = TUrl(url).GetHost();
   fID      = id;
   fMethod  = meth;
   fMethodName = "";
   fOffSet  = offset;
   fToken   = token;
   fUser    = TUrl(url).GetUser();

   // Keep official list updated with active TSecContexts
   if (fOffSet > -1) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSecContexts()->Add(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TSecContext::TSecContext(const TSecContext& sc) :
  TObject(sc),
  fContext(sc.fContext),
  fCleanup(sc.fCleanup),
  fExpDate(sc.fExpDate),
  fHost(sc.fHost),
  fID(sc.fID),
  fMethod(sc.fMethod),
  fMethodName(sc.fMethodName),
  fOffSet(sc.fOffSet),
  fToken(sc.fToken),
  fUser(sc.fUser)
{
}

////////////////////////////////////////////////////////////////////////////////
///assignement operator

TSecContext& TSecContext::operator=(const TSecContext& sc)
{
   if(this!=&sc) {
      TObject::operator=(sc);
      fContext=sc.fContext;
      fCleanup=sc.fCleanup;
      fExpDate=sc.fExpDate;
      fHost=sc.fHost;
      fID=sc.fID;
      fMethod=sc.fMethod;
      fMethodName=sc.fMethodName;
      fOffSet=sc.fOffSet;
      fToken=sc.fToken;
      fUser=sc.fUser;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Dtor: delete (deActivate, local/remote cleanup, list removal)
/// all what is still active

TSecContext::~TSecContext()
{
   Cleanup();
}
////////////////////////////////////////////////////////////////////////////////
/// Cleanup what is still active

void TSecContext::Cleanup()
{
   if (IsActive()) {
      CleanupSecContext(kTRUE);
      DeActivate("R");
      // All have been remotely Deactivated
      TIter nxtl(gROOT->GetListOfSecContexts());
      TSecContext *nscl;
      while ((nscl = (TSecContext *)nxtl())) {
         if (nscl != this && !strcmp(nscl->GetHost(), fHost.Data())) {
            // Need to set ofs=-1 to avoid sending another
            // cleanup request
            nscl->DeActivate("");
         }
      }
   }

   // Delete the cleanup list
   if (fCleanup) {
      fCleanup->Delete();
      delete fCleanup;
      fCleanup = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set OffSet to -1 and expiring Date to default
/// Remove from the list
/// If Opt contains "C" or "c", ask for remote cleanup
/// If Opt contains "R" or "r", remove from the list
/// Default Opt="CR"

void TSecContext::DeActivate(Option_t *Opt)
{
   // Ask remote cleanup of this context
   Bool_t clean = (strstr(Opt,"C") || strstr(Opt,"c"));
   if (clean && fOffSet > -1)
      CleanupSecContext(kFALSE);

   Bool_t remove = (strstr(Opt,"R") || strstr(Opt,"r"));
   if (remove && fOffSet > -1){
      R__LOCKGUARD(gROOTMutex);
      // Remove from the global list
      gROOT->GetListOfSecContexts()->Remove(this);
   }

   // Set inactive
   fOffSet  = -1;
   fExpDate = kROOTTZERO;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new TSecContextCleanup
/// Internally is added to the list

void TSecContext::AddForCleanup(Int_t port, Int_t proto, Int_t type)
{
   TSecContextCleanup *tscc = new TSecContextCleanup(port, proto, type);
   fCleanup->Add(tscc);

}

////////////////////////////////////////////////////////////////////////////////
/// Checks if this security context is for method named 'methname'
/// Case sensitive.

Bool_t TSecContext::IsA(const char *methname)
{
   return Bool_t(!strcmp(methname, GetMethodName()));
}

////////////////////////////////////////////////////////////////////////////////
/// Check remote OffSet and expiring Date

Bool_t TSecContext::IsActive() const
{
   if (fOffSet > -1 && fExpDate > TDatime())
      return kTRUE;
   // Invalid
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// If opt is "F" (default) print object content.
/// If opt is "<number>" print in special form for calls within THostAuth
/// with cardinality <number>
/// If opt is "S" prints short in-line form for calls within TFTP,
/// TSlave, TProof ...

void TSecContext::Print(Option_t *opt) const
{
   char aOrd[16] = {0};
   char aSpc[16] = {0};

   // Check if option is numeric
   Int_t ord = -1, i = 0;
   for (; i < (Int_t)strlen(opt); i++) {
      if (opt[i] < 48 || opt[i] > 57) {
         ord = -2;
         break;
      }
   }
   // If numeric get the cardinality and prepare the strings
   if (ord == -1)
      ord = atoi(opt);

   // If asked to print ordinal number, preapre the string
   if (ord > -1) {
      snprintf(aOrd, sizeof(aOrd), "%d)", ord);
      // and take care of alignment
      Int_t len=strlen(aOrd);
      while (len--)
         strlcat(aSpc, " ", sizeof(aSpc));
   }

   if (!strncasecmp(opt,"F",1)) {
      Info("Print",
           "+------------------------------------------------------+");
      Info("Print",
           "+ Host:%s Method:%d (%s) User:'%s'",
            GetHost(), fMethod, GetMethodName(),
            fUser.Data());
      Info("Print",
           "+         OffSet:%d, id:%s", fOffSet, fID.Data());
      if (fOffSet > -1)
         Info("Print",
           "+         Expiration time: %s",fExpDate.AsString());
      Info("Print",
           "+------------------------------------------------------+");
   } else if (!strncasecmp(opt,"S",1)) {
      if (fOffSet > -1) {
         Printf("Security context:     Method: %d (%s) expiring on %s",
                fMethod, GetMethodName(),
                fExpDate.AsString());
      } else {
         Printf("Security context:     Method: %d (%s) not reusable",
                fMethod, GetMethodName());
      }
   } else {
      // special printing form for THostAuth
      Info("PrintEstblshed","+ %s h:%s met:%d (%s) us:'%s'",
            aOrd, GetHost(), fMethod, GetMethodName(),
            fUser.Data());
      Info("PrintEstblshed","+ %s offset:%d id:%s", aSpc, fOffSet, fID.Data());
      if (fOffSet > -1)
         Info("PrintEstblshed","+ %s expiring: %s",aSpc,fExpDate.AsString());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns short string with relevant information about this
/// security context

const char *TSecContext::AsString(TString &out)
{
   if (fOffSet > -1) {
      char expdate[32];
      out = Form("Method: %d (%s) expiring on %s",
                 fMethod, GetMethodName(), fExpDate.AsString(expdate));
   } else {
      if (fOffSet == -1)
         out = Form("Method: %d (%s) not reusable", fMethod, GetMethodName());
      else if (fOffSet == -3)
         out = Form("Method: %d (%s) authorized by /etc/hosts.equiv or $HOME/.rhosts",
                    fMethod, GetMethodName());
      else if (fOffSet == -4)
         out = Form("No authentication required remotely");
   }
   return out.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Ask remote client to cleanup security context 'ctx'
/// If 'all', all sec context with the same host as ctx
/// are cleaned.

Bool_t TSecContext::CleanupSecContext(Bool_t)
{
   AbstractMethod("CleanupSecContext");
   return kFALSE;
}
