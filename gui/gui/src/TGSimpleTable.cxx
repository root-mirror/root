// Author: Roel Aaij 21/07/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGTable.h"
#include "TGWindow.h"
#include "TGResourcePool.h"
#include "TRandom3.h"
#include "TGSimpleTableInterface.h"
#include "TGSimpleTable.h"

ClassImp(TGSimpleTable);

////////////////////////////////////////////////////////////////////////////////

/* Begin_Html
<center><h2>TGSimpleTable</h2></center>
<br><br>
To provide a simple class to visualize an array of doubles, the class
TGSimpleTable is provided. TGSimpleTable creates it's own
TGSimpleTableInterface. For more information, see the documentation of
TGTable
<br><br>
The interface is accesible through the GetInterface() method.
End_Html
*/

////////////////////////////////////////////////////////////////////////////////
/// TGSimpleTable constuctor.

TGSimpleTable::TGSimpleTable(TGWindow *p, Int_t id, Double_t **data,
                             UInt_t nrows, UInt_t ncolumns)
   : TGTable(p, id, 0, nrows, ncolumns)
{
   TGSimpleTableInterface *iface = new TGSimpleTableInterface(data, nrows,
                                                              ncolumns);
   SetInterface(iface,nrows, ncolumns);
}

////////////////////////////////////////////////////////////////////////////////
/// TGSimpleTable destructor.

TGSimpleTable::~TGSimpleTable()
{
   delete fInterface;
}

