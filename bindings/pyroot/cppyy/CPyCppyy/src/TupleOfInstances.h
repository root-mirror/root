#ifndef CPYCPPYY_TUPLEOFINSTANCES_H
#define CPYCPPYY_TUPLEOFINSTANCES_H

namespace CPyCppyy {

/** Representation of C-style array of instances
      @author  WLAV
      @date    02/10/2014
      @version 1.0
 */

//- custom tuple type that can pass through C-style arrays -------------------
extern PyTypeObject TupleOfInstances_Type;
extern PyTypeObject InstanceArrayIter_Type;

template<typename T>
inline bool TupleOfInstances_Check(T* object)
{
    return object && PyObject_TypeCheck(object, &TupleOfInstances_Type);
}

template<typename T>
inline bool TupleOfInstances_CheckExact(T* object)
{
    return object && Py_TYPE(object) == &TupleOfInstances_Type;
}

PyObject* TupleOfInstances_New(
    Cppyy::TCppObject_t address, Cppyy::TCppType_t klass, dim_t ndims, dims_t dims);

} // namespace CPyCppyy

#endif // !CPYCPPYY_TUPLEOFINSTANCES_H
