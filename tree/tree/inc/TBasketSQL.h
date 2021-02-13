// @(#)root/tree:$Id$
// Author: Philippe Canal 2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef TBASKETSQL_H
#define TBASKETSQL_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBasketSQL                                                           //
//                                                                      //
// Implement TBasket for a SQL backend                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include "TBasket.h"

class TSQLResult;
class TSQLRow;
class TBufferSQL;

class TBasketSQL : public TBasket
{

private:
   TBasketSQL(const TBasketSQL&);            // TBasketSQL objects are not copiable.
   TBasketSQL& operator=(const TBasketSQL&); // TBasketSQL objects are not copiable.

protected:
   TSQLResult **fResultPtr;    //!
   TSQLRow    **fRowPtr;       //!
   TString      *fInsertQuery; //!

public:
   TBasketSQL();
   TBasketSQL(const char *name, const char *title,
              TBranch *branch, TSQLResult **rs,
              TString *insert_query, std::vector<Int_t> *vc, TSQLRow **row);
   ~TBasketSQL();
   void    PrepareBasket(Long64_t entry);
   virtual Int_t   ReadBasketBuffers(Long64_t pos, Int_t len, TFile *file);
   virtual Int_t   ReadBasketBytes(Long64_t pos, TFile *file);

   TSQLResult * GetResultSet() { return *fResultPtr;}
   void CreateBuffer(const char *name, TString title, std::vector<Int_t> * vc, TBranch *branch, TSQLResult ** rs);

   void Update(Int_t offset, Int_t skipped);

   ClassDef(TBasketSQL,1)  //the TBranch buffers

};

#endif
