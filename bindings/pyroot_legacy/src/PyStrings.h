// Author: Wim Lavrijsen, Nov 2008

#ifndef PYROOT_PYSTRINGS_H
#define PYROOT_PYSTRINGS_H

// ROOT
#include "DllImport.h"


namespace PyROOT {

// python strings kept for performance reasons

   namespace PyStrings {

      R__EXTERN PyObject* gBases;
      R__EXTERN PyObject* gBase;
      R__EXTERN PyObject* gClass;
      R__EXTERN PyObject* gCppEq;
      R__EXTERN PyObject* gCppNe;
      R__EXTERN PyObject* gDeref;
      R__EXTERN PyObject* gDict;
      R__EXTERN PyObject* gEmptyString;
      R__EXTERN PyObject* gEq;
      R__EXTERN PyObject* gFollow;
      R__EXTERN PyObject* gGetItem;
      R__EXTERN PyObject* gInit;
      R__EXTERN PyObject* gIter;
      R__EXTERN PyObject* gLen;
      R__EXTERN PyObject* gLifeLine;
      R__EXTERN PyObject* gModule;
      R__EXTERN PyObject* gMRO;
      R__EXTERN PyObject* gName;
      R__EXTERN PyObject* gCppName;
      R__EXTERN PyObject* gCppNameNew;
      R__EXTERN PyObject* gNe;
      R__EXTERN PyObject* gTypeCode;

      R__EXTERN PyObject* gAdd;
      R__EXTERN PyObject* gSub;
      R__EXTERN PyObject* gMul;
      R__EXTERN PyObject* gDiv;

      R__EXTERN PyObject* gAt;
      R__EXTERN PyObject* gBegin;
      R__EXTERN PyObject* gEnd;
      R__EXTERN PyObject* gFirst;
      R__EXTERN PyObject* gSecond;
      R__EXTERN PyObject* gSize;
      R__EXTERN PyObject* gGetSize;
      R__EXTERN PyObject* ggetSize;
      R__EXTERN PyObject* gTemplate;
      R__EXTERN PyObject* gVectorAt;

      R__EXTERN PyObject* gBranch;
      R__EXTERN PyObject* gFitFCN;
      R__EXTERN PyObject* gROOTns;
      R__EXTERN PyObject* gSetBranchAddress;
      R__EXTERN PyObject* gSetFCN;
      R__EXTERN PyObject* gTClassDynCast;

   } // namespace PyStrings

   Bool_t CreatePyStrings();
   PyObject* DestroyPyStrings();

} // namespace PyROOT

#endif // !PYROOT_PYSTRINGS_H
