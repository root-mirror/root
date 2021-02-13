// Author: Roel Aaij 30/08/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGTreeTable.h"
#include "TTreeTableInterface.h"

ClassImp(TGTreeTable);

/** \class TGTreeTable

TGTreeTable is a TGTable that owns it's own interface.
It can be used to view a TTree. If an expression is given to the
constructor, it will be used to define the columns. A selection can
also be given. This selection is applied to the TTree as a
TEntryList. See the documentation of TGTable for more information

The interface is accesible after the creation through the
GetInterface() method.
*/

////////////////////////////////////////////////////////////////////////////////
/// TGTreeTable constructor.

TGTreeTable::TGTreeTable(TGWindow *p, Int_t id, TTree *tree,
                         const char *expression, const char *selection,
                         const char *option, UInt_t nrows, UInt_t ncolumns)
   : TGTable(p, id, 0, nrows, ncolumns)
{
   TTreeTableInterface *iface = new TTreeTableInterface(tree, expression,
                                                        selection, option);
   SetInterface(iface, nrows, ncolumns);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// TGTreeTable destructor.

TGTreeTable::~TGTreeTable()
{
   //FIXME this causes a double delete segfault, why???
//    delete fInterface;
}

