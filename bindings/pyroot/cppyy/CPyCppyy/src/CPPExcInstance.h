#ifndef CPYCPPYY_CPPEXCINSTANCE_H
#define CPYCPPYY_CPPEXCINSTANCE_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// CpyCppyy::CPPExceptionInstance                                           //
//                                                                          //
// Python-side proxy, encapsulaties a C++ exception object.                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace CPyCppyy {

class CPPExcInstance {
public:
    PyBaseExceptionObject fBase;
    PyObject*             fCppInstance;
    PyObject*             fTopMessage;
};


//- object proxy type and type verification ----------------------------------
CPYCPPYY_IMPORT PyTypeObject CPPExcInstance_Type;

template<typename T>
inline bool CPPExcInstance_Check(T* object)
{
    return object && PyObject_TypeCheck(object, &CPPExcInstance_Type);
}

template<typename T>
inline bool CPPExcInstance_CheckExact(T* object)
{
    return object && Py_TYPE(object) == &CPPExcInstance_Type;
}

} // namespace CPyCppyy

#endif // !CPYCPPYY_CPPEXCINSTANCE_H
