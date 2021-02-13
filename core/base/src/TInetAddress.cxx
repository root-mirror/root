// @(#)root/base:$Id$
// Author: Fons Rademakers   16/12/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TInetAddress
\ingroup Base

This class represents an Internet Protocol (IP) address.
*/

#include "TInetAddress.h"
#include "TBuffer.h"

ClassImp(TInetAddress);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor. Used in case of unknown host. Not a valid address.

TInetAddress::TInetAddress()
{
   fHostname  = "UnknownHost";
   AddAddress(0);
   fFamily    = -1;
   fPort      = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Create TInetAddress. Private ctor. TInetAddress objects can only
/// be created via the friend classes TSystem, TServerSocket and TSocket.
/// Use the IsValid() method to check the validity of a TInetAddress.

TInetAddress::TInetAddress(const char *host, UInt_t addr, Int_t family, Int_t port)
{
   AddAddress(addr);
   if (!strcmp(host, "????") || !strcmp(host, "UnNamedHost"))
      fHostname = GetHostAddress();
   else
      fHostname = host;
   fFamily    = family;
   fPort      = port;
}

////////////////////////////////////////////////////////////////////////////////
/// TInetAddress copy ctor.

TInetAddress::TInetAddress(const TInetAddress &adr) : TObject(adr)
{
   fHostname  = adr.fHostname;
   fFamily    = adr.fFamily;
   fPort      = adr.fPort;
   fAddresses = adr.fAddresses;
   fAliases   = adr.fAliases;
}

////////////////////////////////////////////////////////////////////////////////
/// TInetAddress assignment operator.

TInetAddress& TInetAddress::operator=(const TInetAddress &rhs)
{
   if (this != &rhs) {
      TObject::operator=(rhs);
      fHostname  = rhs.fHostname;
      fFamily    = rhs.fFamily;
      fPort      = rhs.fPort;
      fAddresses = rhs.fAddresses;
      fAliases   = rhs.fAliases;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the raw IP address in host byte order. The highest
/// order byte position is in addr[0]. To be prepared for 64-bit
/// IP addresses an array of bytes is returned.
/// User must delete allocated memory.

UChar_t *TInetAddress::GetAddressBytes() const
{
   UChar_t *addr = new UChar_t[4];

   addr[0] = (UChar_t) ((fAddresses[0] >> 24) & 0xFF);
   addr[1] = (UChar_t) ((fAddresses[0] >> 16) & 0xFF);
   addr[2] = (UChar_t) ((fAddresses[0] >> 8) & 0xFF);
   addr[3] = (UChar_t)  (fAddresses[0] & 0xFF);

   return addr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the IP address string "%d.%d.%d.%d", use it to convert
/// alternative addresses obtained via GetAddresses().
/// Copy string immediately, it will be reused. Static function.

const char *TInetAddress::GetHostAddress(UInt_t addr)
{
   return Form("%d.%d.%d.%d", (addr >> 24) & 0xFF,
                              (addr >> 16) & 0xFF,
                              (addr >>  8) & 0xFF,
                               addr & 0xFF);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the IP address string "%d.%d.%d.%d".
/// Copy string immediately, it will be reused.

const char *TInetAddress::GetHostAddress() const
{
   return GetHostAddress(fAddresses[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Print internet address as string.

void TInetAddress::Print(Option_t *) const
{
   if (fPort == -1)
      Printf("%s/%s (not connected)", GetHostName(), GetHostAddress());
   else
      Printf("%s/%s (port %d)", GetHostName(), GetHostAddress(), fPort);

   int i = 0;
   AddressList_t::const_iterator ai;
   for (ai = fAddresses.begin(); ai != fAddresses.end(); ++ai) {
      if (!i) printf("%s:", fAddresses.size() == 1 ? "Address" : "Addresses");
      printf(" %s", GetHostAddress(*ai));
      i++;
   }
   if (i) printf("\n");

   i = 0;
   AliasList_t::const_iterator ali;
   for (ali = fAliases.begin(); ali != fAliases.end(); ++ali) {
      if (!i) printf("%s:", fAliases.size() == 1 ? "Alias" : "Aliases");
      printf(" %s", ali->Data());
      i++;
   }
   if (i) printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Add alternative address to list of addresses.

void TInetAddress::AddAddress(UInt_t addr)
{
   fAddresses.push_back(addr);
}

////////////////////////////////////////////////////////////////////////////////
/// Add alias to list of aliases.

void TInetAddress::AddAlias(const char *alias)
{
   fAliases.push_back(TString(alias));
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TInetAddress.

void TInetAddress::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 2) {
         R__b.ReadClassBuffer(TInetAddress::Class(), this, R__v, R__s, R__c);
         return;
      }
      // process old versions before automatic schema evolution
      UInt_t address;
      TObject::Streamer(R__b);
      fHostname.Streamer(R__b);
      R__b >> address;
      R__b >> fFamily;
      R__b >> fPort;
      if (R__v == 1)
         fAddresses.push_back(address);
      if (R__v > 1) {
         TInetAddress::AddressList_t &R__stl1 =  fAddresses;
         R__stl1.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl1.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            unsigned int R__t1;
            R__b >> R__t1;
            R__stl1.push_back(R__t1);
         }
         TInetAddress::AliasList_t &R__stl2 =  fAliases;
         R__stl2.clear();
         R__b >> R__n;
         R__stl2.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            TString R__t2;
            R__t2.Streamer(R__b);
            R__stl2.push_back(R__t2);
         }
      }
      R__b.CheckByteCount(R__s, R__c, TInetAddress::IsA());
   } else {
      R__b.WriteClassBuffer(TInetAddress::Class(), this);
   }
}
