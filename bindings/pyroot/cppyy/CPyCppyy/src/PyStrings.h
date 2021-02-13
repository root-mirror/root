#ifndef CPYCPPYY_PYSTRINGS_H
#define CPYCPPYY_PYSTRINGS_H

namespace CPyCppyy {

// python strings kept for performance reasons

namespace PyStrings {

    extern PyObject* gAssign;
    extern PyObject* gBases;
    extern PyObject* gBase;
    extern PyObject* gCppName;
    extern PyObject* gDeref;
    extern PyObject* gPreInc;
    extern PyObject* gPostInc;
    extern PyObject* gDict;
    extern PyObject* gEmptyString;
    extern PyObject* gEq;
    extern PyObject* gFollow;
    extern PyObject* gGetItem;
    extern PyObject* gGetNoCheck;
    extern PyObject* gInit;
    extern PyObject* gIter;
    extern PyObject* gLen;
    extern PyObject* gLifeLine;
    extern PyObject* gModule;
    extern PyObject* gMRO;
    extern PyObject* gName;
    extern PyObject* gNe;
    extern PyObject* gTypeCode;
    extern PyObject* gCTypesType;

    extern PyObject* gUnderlying;

    extern PyObject* gAdd;
    extern PyObject* gSub;
    extern PyObject* gMul;
    extern PyObject* gDiv;

    extern PyObject* gLShift;
    extern PyObject* gLShiftC;

    extern PyObject* gAt;
    extern PyObject* gBegin;
    extern PyObject* gEnd;
    extern PyObject* gFirst;
    extern PyObject* gSecond;
    extern PyObject* gSize;
    extern PyObject* gTemplate;
    extern PyObject* gVectorAt;

    extern PyObject* gCppReal;
    extern PyObject* gCppImag;

    extern PyObject* gThisModule;

    extern PyObject* gNoImplicit;
    extern PyObject* gDispInit;

    extern PyObject* gExPythonize;
    extern PyObject* gPythonize;

} // namespace PyStrings

bool CreatePyStrings();
PyObject* DestroyPyStrings();

} // namespace CPyCppyy

#endif // !CPYCPPYY_PYSTRINGS_H
