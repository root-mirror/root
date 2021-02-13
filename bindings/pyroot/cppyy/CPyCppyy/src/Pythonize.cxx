// Bindings
#include "CPyCppyy.h"
#include "Pythonize.h"
#include "Converters.h"
#include "CPPInstance.h"
#include "CPPOverload.h"
#include "CustomPyTypes.h"
#include "LowLevelViews.h"
#include "ProxyWrappers.h"
#include "PyCallable.h"
#include "PyStrings.h"
#include "TypeManip.h"
#include "Utility.h"

// Standard
#include <algorithm>
#include <complex>
#include <set>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>


//- data and local helpers ---------------------------------------------------
namespace CPyCppyy {
    extern PyObject* gThisModule;
    extern std::map<std::string, std::vector<PyObject*>> gPythonizations;
}

namespace {

// for convenience
using namespace CPyCppyy;

//-----------------------------------------------------------------------------
bool HasAttrDirect(PyObject* pyclass, PyObject* pyname, bool mustBeCPyCppyy = false) {
// prevents calls to Py_TYPE(pyclass)->tp_getattr, which is unnecessary for our
// purposes here and could tickle problems w/ spurious lookups into ROOT meta
    PyObject* dct = PyObject_GetAttr(pyclass, PyStrings::gDict);
    if (dct) {
        PyObject* attr = PyObject_GetItem(dct, pyname);
        Py_DECREF(dct);
        if (attr) {
            bool ret = !mustBeCPyCppyy || CPPOverload_Check(attr);
            Py_DECREF(attr);
            return ret;
        }
    }
    PyErr_Clear();
    return false;
}

//-----------------------------------------------------------------------------
inline bool IsTemplatedSTLClass(const std::string& name, const std::string& klass) {
// Scan the name of the class and determine whether it is a template instantiation.
    auto pos = name.find(klass);
    return (pos == 0 || pos == 5) && name.find("::", name.rfind(">")) == std::string::npos;
}

// to prevent compiler warnings about const char* -> char*
inline PyObject* CallPyObjMethod(PyObject* obj, const char* meth)
{
// Helper; call method with signature: obj->meth().
    Py_INCREF(obj);
    PyObject* result = PyObject_CallMethod(obj, const_cast<char*>(meth), const_cast<char*>(""));
    Py_DECREF(obj);
    return result;
}

//-----------------------------------------------------------------------------
inline PyObject* CallPyObjMethod(PyObject* obj, const char* meth, PyObject* arg1)
{
// Helper; call method with signature: obj->meth(arg1).
    Py_INCREF(obj);
    PyObject* result = PyObject_CallMethod(
        obj, const_cast<char*>(meth), const_cast<char*>("O"), arg1);
    Py_DECREF(obj);
    return result;
}

//-----------------------------------------------------------------------------
PyObject* PyStyleIndex(PyObject* self, PyObject* index)
{
// Helper; converts python index into straight C index.
    Py_ssize_t idx = PyInt_AsSsize_t(index);
    if (idx == (Py_ssize_t)-1 && PyErr_Occurred())
        return nullptr;

    Py_ssize_t size = PySequence_Size(self);
    if (idx >= size || (idx < 0 && idx < -size)) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return nullptr;
    }

    PyObject* pyindex = nullptr;
    if (idx >= 0) {
        Py_INCREF(index);
        pyindex = index;
    } else
        pyindex = PyLong_FromSsize_t(size+idx);

    return pyindex;
}

//-----------------------------------------------------------------------------
inline bool AdjustSlice(const Py_ssize_t nlen, Py_ssize_t& start, Py_ssize_t& stop, Py_ssize_t& step)
{
// Helper; modify slice range to match the container.
    if ((step > 0 && stop <= start) || (step < 0 && start <= stop))
        return false;

    if (start < 0) start = 0;
    if (start >= nlen) start = nlen-1;
    if (step >= nlen) step = nlen;

    stop = step > 0 ? std::min(nlen, stop) : (stop >= 0 ? stop : -1);
    return true;
}

//-----------------------------------------------------------------------------
inline PyObject* CallSelfIndex(CPPInstance* self, PyObject* idx, PyObject* pymeth)
{
// Helper; call method with signature: meth(pyindex).
    Py_INCREF((PyObject*)self);
    PyObject* pyindex = PyStyleIndex((PyObject*)self, idx);
    if (!pyindex) {
        Py_DECREF((PyObject*)self);
        return nullptr;
    }

    PyObject* result = PyObject_CallMethodObjArgs((PyObject*)self, pymeth, pyindex, nullptr);
    Py_DECREF(pyindex);
    Py_DECREF((PyObject*)self);
    return result;
}

//- "smart pointer" behavior ---------------------------------------------------
PyObject* DeRefGetAttr(PyObject* self, PyObject* name)
{
// Follow operator*() if present (available in python as __deref__), so that
// smart pointers behave as expected.
    if (name == PyStrings::gTypeCode || name == PyStrings::gCTypesType) {
    // TODO: these calls come from TemplateProxy and are unlikely to be needed in practice,
    // whereas as-is, they can accidentally dereference the result of end() on some STL
    // containers. Obviously, this is a dumb hack that should be resolved more fundamentally.
        PyErr_SetString(PyExc_AttributeError, CPyCppyy_PyText_AsString(name));
        return nullptr;
    }

    if (!CPyCppyy_PyText_Check(name))
        PyErr_SetString(PyExc_TypeError, "getattr(): attribute name must be string");

    PyObject* pyptr = PyObject_CallMethodObjArgs(self, PyStrings::gDeref, nullptr);
    if (!pyptr)
        return nullptr;

// prevent a potential infinite loop
    if (Py_TYPE(pyptr) == Py_TYPE(self)) {
        PyObject* val1 = PyObject_Str(self);
        PyObject* val2 = PyObject_Str(name);
        PyErr_Format(PyExc_AttributeError, "%s has no attribute \'%s\'",
            CPyCppyy_PyText_AsString(val1), CPyCppyy_PyText_AsString(val2));
        Py_DECREF(val2);
        Py_DECREF(val1);

        Py_DECREF(pyptr);
        return nullptr;
    }

    PyObject* result = PyObject_GetAttr(pyptr, name);
    Py_DECREF(pyptr);
    return result;
}

//-----------------------------------------------------------------------------
PyObject* FollowGetAttr(PyObject* self, PyObject* name)
{
// Follow operator->() if present (available in python as __follow__), so that
// smart pointers behave as expected.
    if (!CPyCppyy_PyText_Check(name))
        PyErr_SetString(PyExc_TypeError, "getattr(): attribute name must be string");

    PyObject* pyptr = PyObject_CallMethodObjArgs(self, PyStrings::gFollow, nullptr);
    if (!pyptr)
         return nullptr;

    PyObject* result = PyObject_GetAttr(pyptr, name);
    Py_DECREF(pyptr);
    return result;
}


//- vector behavior as primitives ----------------------------------------------
#if PY_VERSION_HEX < 0x03040000
#define PyObject_LengthHint _PyObject_LengthHint
#endif

// TODO: can probably use the below getters in the InitializerListConverter
struct ItemGetter {
    ItemGetter(PyObject* pyobj) : fPyObject(pyobj) { Py_INCREF(fPyObject); }
    virtual ~ItemGetter() { Py_DECREF(fPyObject); }
    virtual Py_ssize_t size() = 0;
    virtual PyObject* get() = 0;
    PyObject* fPyObject;
};

struct CountedItemGetter : public ItemGetter {
    CountedItemGetter(PyObject* pyobj) : ItemGetter(pyobj), fCur(0) {}
    Py_ssize_t fCur;
};

struct TupleItemGetter : public CountedItemGetter {
    using CountedItemGetter::CountedItemGetter;
    virtual Py_ssize_t size() { return PyTuple_GET_SIZE(fPyObject); }
    virtual PyObject* get() {
        if (fCur < PyTuple_GET_SIZE(fPyObject)) {
            PyObject* item = PyTuple_GET_ITEM(fPyObject, fCur++);
            Py_INCREF(item);
            return item;
        }
        PyErr_SetString(PyExc_StopIteration, "end of tuple");
        return nullptr;
    }
};

struct ListItemGetter : public CountedItemGetter {
    using CountedItemGetter::CountedItemGetter;
    virtual Py_ssize_t size() { return PyList_GET_SIZE(fPyObject); }
    virtual PyObject* get() {
        if (fCur < PyList_GET_SIZE(fPyObject)) {
            PyObject* item = PyList_GET_ITEM(fPyObject, fCur++);
            Py_INCREF(item);
            return item;
        }
        PyErr_SetString(PyExc_StopIteration, "end of list");
        return nullptr;
    }
};

struct SequenceItemGetter : public CountedItemGetter {
    using CountedItemGetter::CountedItemGetter;
    virtual Py_ssize_t size() {
        Py_ssize_t sz = PySequence_Size(fPyObject);
        if (sz < 0) {
            PyErr_Clear();
            return PyObject_LengthHint(fPyObject, 8);
        }
        return sz;
    }
    virtual PyObject* get() { return PySequence_GetItem(fPyObject, fCur++); }
};

struct IterItemGetter : public ItemGetter {
    using ItemGetter::ItemGetter;
    virtual Py_ssize_t size() { return PyObject_LengthHint(fPyObject, 8); }
    virtual PyObject* get() { return (*(Py_TYPE(fPyObject)->tp_iternext))(fPyObject); }
};

PyObject* VectorInit(PyObject* self, PyObject* args, PyObject* /* kwds */)
{
// Specialized vector constructor to allow construction from containers; allowing
// such construction from initializer_list instead would possible, but can be
// error-prone. This use case is common enough for std::vector to implement it
// directly, except for arrays (which can be passed wholesale) and strings (which
// won't convert properly as they'll be seen as buffers)

    ItemGetter* getter = nullptr;
    if (PyTuple_GET_SIZE(args) == 1) {
        PyObject* fi = PyTuple_GET_ITEM(args, 0);
        if (CPyCppyy_PyText_Check(fi) || PyBytes_Check(fi)) {
            PyErr_SetString(PyExc_TypeError, "can not convert string to vector");
            return nullptr;
        }
    // TODO: this only tests for new-style buffers, which is too strict, but a
    // generic check for Py_TYPE(fi)->tp_as_buffer is too loose (note that the
    // main use case is numpy, which offers the new interface)
        if (!PyObject_CheckBuffer(fi)) {
            if (PyTuple_CheckExact(fi))
                getter = new TupleItemGetter(fi);
            else if (PyList_CheckExact(fi))
                getter = new ListItemGetter(fi);
            else if (PySequence_Check(fi))
                getter = new SequenceItemGetter(fi);
            else {
                PyObject* iter = PyObject_GetIter(fi);
                if (iter) {
                    getter = new IterItemGetter{iter};
                    Py_DECREF(iter);
                }
                else PyErr_Clear();
            }
        }
    }

    if (getter) {
    // construct an empty vector, then back-fill it
        PyObject* mname = CPyCppyy_PyText_FromString("__real_init");
        PyObject* result = PyObject_CallMethodObjArgs(self, mname, nullptr);
        Py_DECREF(mname);
        if (!result) {
            delete getter;
            return result;
        }

        Py_ssize_t sz = getter->size();
        if (sz < 0) {
            delete getter;
            return nullptr;
        }

    // reserve memory as appliable
        if (0 < sz) {
            PyObject* res = PyObject_CallMethod(self, (char*)"reserve", (char*)"n", sz);
            Py_DECREF(res);
        } else { // empty container
            delete getter;
            return result;
        }

        bool fill_ok = true;

    // two main options: a list of lists (or tuples), or a list of objects; the former
    // are emplace_back'ed, the latter push_back'ed
        PyObject* fi = PySequence_GetItem(PyTuple_GET_ITEM(args, 0), 0);
        if (!fi) PyErr_Clear();
        if (fi && (PyTuple_CheckExact(fi) || PyList_CheckExact(fi))) {
        // use emplace_back to construct the vector entries one by one
            PyObject* eb_call = PyObject_GetAttrString(self, (char*)"emplace_back");
            PyObject* vtype = PyObject_GetAttrString((PyObject*)Py_TYPE(self), "value_type");
            bool value_is_vector = false;
            if (vtype && CPyCppyy_PyText_Check(vtype)) {
            // if the value_type is a vector, then allow for initialization from sequences
                if (std::string(CPyCppyy_PyText_AsString(vtype)).rfind("std::vector", 0) != std::string::npos)
                    value_is_vector = true;
            } else
                PyErr_Clear();
            Py_XDECREF(vtype);

            if (eb_call) {
                PyObject* eb_args;
                for (int i = 0; /* until break */; ++i) {
                    PyObject* item = getter->get();
                    if (item) {
                        if (value_is_vector && PySequence_Check(item)) {
                            eb_args = PyTuple_New(1);
                            PyTuple_SET_ITEM(eb_args, 0, item);
                        } else if (PyTuple_CheckExact(item)) {
                            eb_args = item;
                        } else if (PyList_CheckExact(item)) {
                            Py_ssize_t isz = PyList_GET_SIZE(item);
                            eb_args = PyTuple_New(isz);
                            for (Py_ssize_t j = 0; j < isz; ++j) {
                                PyObject* iarg = PyList_GET_ITEM(item, j);
                                Py_INCREF(iarg);
                                PyTuple_SET_ITEM(eb_args, j, iarg);
                            }
                            Py_DECREF(item);
                        } else {
                            Py_DECREF(item);
                            PyErr_Format(PyExc_TypeError, "argument %d is not a tuple or list", i);
                            fill_ok = false;
                            break;
                        }
                        PyObject* ebres = PyObject_CallObject(eb_call, eb_args);
                        Py_DECREF(eb_args);
                        if (!ebres) {
                            fill_ok = false;
                            break;
                        }
                        Py_DECREF(ebres);
                    } else {
                        if (PyErr_Occurred()) {
                            if (!(PyErr_ExceptionMatches(PyExc_IndexError) ||
                                  PyErr_ExceptionMatches(PyExc_StopIteration)))
                                fill_ok = false;
                            else { PyErr_Clear(); }
                        }
                        break;
                    }
                }
                Py_DECREF(eb_call);
            }
        } else {
        // use push_back to add the vector entries one by one
            PyObject* pb_call = PyObject_GetAttrString(self, (char*)"push_back");
            if (pb_call) {
                for (;;) {
                    PyObject* item = getter->get();
                    if (item) {
                        PyObject* pbres = PyObject_CallFunctionObjArgs(pb_call, item, nullptr);
                        Py_DECREF(item);
                        if (!pbres) {
                            fill_ok = false;
                            break;
                        }
                        Py_DECREF(pbres);
                    } else {
                        if (PyErr_Occurred()) {
                            if (!(PyErr_ExceptionMatches(PyExc_IndexError) ||
                                  PyErr_ExceptionMatches(PyExc_StopIteration)))
                                fill_ok = false;
                            else { PyErr_Clear(); }
                        }
                        break;
                    }
                }
                Py_DECREF(pb_call);
            }
        }
        Py_XDECREF(fi);
        delete getter;

        if (!fill_ok) {
            Py_DECREF(result);
            return nullptr;
        }

        return result;
    }

// The given argument wasn't iterable: simply forward to regular constructor
    PyObject* realInit = PyObject_GetAttrString(self, "__real_init");
    if (realInit) {
        PyObject* result = PyObject_Call(realInit, args, nullptr);
        Py_DECREF(realInit);
        return result;
    }

    return nullptr;
}

//---------------------------------------------------------------------------
PyObject* VectorData(PyObject* self, PyObject*)
{
    PyObject* pydata = CallPyObjMethod(self, "__real_data");
    if (!LowLevelView_Check(pydata)) return pydata;

    PyObject* pylen = PyObject_CallMethodObjArgs(self, PyStrings::gSize, nullptr);
    if (!pylen) {
        PyErr_Clear();
        return pydata;
    }

    long clen = PyInt_AsLong(pylen);
    Py_DECREF(pylen);

// TODO: should be a LowLevelView helper
    Py_buffer& bi = ((LowLevelView*)pydata)->fBufInfo;
    bi.len = clen * bi.itemsize;
    if (bi.ndim == 1 && bi.shape)
        bi.shape[0] = clen;

    return pydata;
}


//-----------------------------------------------------------------------------
static PyObject* vector_iter(PyObject* v) {
    vectoriterobject* vi = PyObject_GC_New(vectoriterobject, &VectorIter_Type);
    if (!vi) return nullptr;

    Py_INCREF(v);
    vi->ii_container = v;
    vi->vi_flags = v->ob_refcnt <= 2 ? 1 : 0;    // 2, b/c of preceding INCREF

    PyObject* pyvalue_type = PyObject_GetAttrString((PyObject*)Py_TYPE(v), "value_type");
    PyObject* pyvalue_size = PyObject_GetAttrString((PyObject*)Py_TYPE(v), "value_size");

    vi->vi_klass = 0;
    if (pyvalue_type && pyvalue_size) {
        PyObject* pydata = CallPyObjMethod(v, "data");
        if (!pydata || Utility::GetBuffer(pydata, '*', 1, vi->vi_data, false) == 0) {
            if (CPPInstance_Check(pydata)) {
                vi->vi_data = ((CPPInstance*)pydata)->GetObjectRaw();
                vi->vi_klass = ((CPPInstance*)pydata)->ObjectIsA(false);
            } else
                vi->vi_data = nullptr;
        }
        Py_XDECREF(pydata);

        vi->vi_converter = vi->vi_klass ? nullptr : CPyCppyy::CreateConverter(CPyCppyy_PyText_AsString(pyvalue_type));
        vi->vi_stride    = PyLong_AsLong(pyvalue_size);
    } else {
        PyErr_Clear();
        vi->vi_data      = nullptr;
        vi->vi_converter = nullptr;
        vi->vi_stride    = 0;
    }

    Py_XDECREF(pyvalue_size);
    Py_XDECREF(pyvalue_type);

    vi->ii_pos = 0;
    vi->ii_len = PySequence_Size(v);

    PyObject_GC_Track(vi);
    return (PyObject*)vi;
}

PyObject* VectorGetItem(CPPInstance* self, PySliceObject* index)
{
// Implement python's __getitem__ for std::vector<>s.
    if (PySlice_Check(index)) {
        if (!self->GetObject()) {
            PyErr_SetString(PyExc_TypeError, "unsubscriptable object");
            return nullptr;
        }

        PyObject* pyclass = (PyObject*)Py_TYPE((PyObject*)self);
        PyObject* nseq = PyObject_CallObject(pyclass, nullptr);

        Py_ssize_t start, stop, step;
        PySlice_GetIndices((CPyCppyy_PySliceCast)index, PyObject_Length((PyObject*)self), &start, &stop, &step);

        const Py_ssize_t nlen = PySequence_Size((PyObject*)self);
        if (!AdjustSlice(nlen, start, stop, step))
            return nseq;

        const Py_ssize_t sign = step < 0 ? -1 : 1;
        for (Py_ssize_t i = start; i*sign < stop*sign; i += step) {
            PyObject* pyidx = PyInt_FromSsize_t(i);
            PyObject* item = PyObject_CallMethodObjArgs((PyObject*)self, PyStrings::gGetNoCheck, pyidx, nullptr);
            CallPyObjMethod(nseq, "push_back", item);
            Py_DECREF(item);
            Py_DECREF(pyidx);
        }

        return nseq;
    }

    return CallSelfIndex(self, (PyObject*)index, PyStrings::gGetNoCheck);
}


static Cppyy::TCppType_t sVectorBoolTypeID = (Cppyy::TCppType_t)0;

PyObject* VectorBoolGetItem(CPPInstance* self, PyObject* idx)
{
// std::vector<bool> is a special-case in C++, and its return type depends on
// the compiler: treat it special here as well
    if (!CPPInstance_Check(self) || self->ObjectIsA() != sVectorBoolTypeID) {
        PyErr_Format(PyExc_TypeError,
            "require object of type std::vector<bool>, but %s given",
            Cppyy::GetScopedFinalName(self->ObjectIsA()).c_str());
        return nullptr;
    }

    if (!self->GetObject()) {
        PyErr_SetString(PyExc_TypeError, "unsubscriptable object");
        return nullptr;
    }

    if (PySlice_Check(idx)) {
        PyObject* pyclass = (PyObject*)Py_TYPE((PyObject*)self);
        PyObject* nseq = PyObject_CallObject(pyclass, nullptr);

        Py_ssize_t start, stop, step;
        PySlice_GetIndices((CPyCppyy_PySliceCast)idx, PyObject_Length((PyObject*)self), &start, &stop, &step);
        const Py_ssize_t nlen = PySequence_Size((PyObject*)self);
        if (!AdjustSlice(nlen, start, stop, step))
            return nseq;

        const Py_ssize_t sign = step < 0 ? -1 : 1;
        for (Py_ssize_t i = start; i*sign < stop*sign; i += step) {
            PyObject* pyidx = PyInt_FromSsize_t(i);
            PyObject* item = PyObject_CallMethodObjArgs((PyObject*)self, PyStrings::gGetItem, pyidx, nullptr);
            CallPyObjMethod(nseq, "push_back", item);
            Py_DECREF(item);
            Py_DECREF(pyidx);
        }

        return nseq;
    }

    PyObject* pyindex = PyStyleIndex((PyObject*)self, idx);
    if (!pyindex)
        return nullptr;

    int index = (int)PyLong_AsLong(pyindex);
    Py_DECREF(pyindex);

// get hold of the actual std::vector<bool> (no cast, as vector is never a base)
    std::vector<bool>* vb = (std::vector<bool>*)self->GetObject();

// finally, return the value
    if (bool((*vb)[index]))
        Py_RETURN_TRUE;
    Py_RETURN_FALSE;
}

PyObject* VectorBoolSetItem(CPPInstance* self, PyObject* args)
{
// std::vector<bool> is a special-case in C++, and its return type depends on
// the compiler: treat it special here as well
    if (!CPPInstance_Check(self) || self->ObjectIsA() != sVectorBoolTypeID) {
        PyErr_Format(PyExc_TypeError,
            "require object of type std::vector<bool>, but %s given",
            Cppyy::GetScopedFinalName(self->ObjectIsA()).c_str());
        return nullptr;
    }

    if (!self->GetObject()) {
        PyErr_SetString(PyExc_TypeError, "unsubscriptable object");
        return nullptr;
    }

    int bval = 0; PyObject* idx = nullptr;
    if (!PyArg_ParseTuple(args, const_cast<char*>("Oi:__setitem__"), &idx, &bval))
        return nullptr;

    PyObject* pyindex = PyStyleIndex((PyObject*)self, idx);
    if (!pyindex)
        return nullptr;

    int index = (int)PyLong_AsLong(pyindex);
    Py_DECREF(pyindex);

// get hold of the actual std::vector<bool> (no cast, as vector is never a base)
    std::vector<bool>* vb = (std::vector<bool>*)self->GetObject();

// finally, set the value
    (*vb)[index] = (bool)bval;

    Py_RETURN_NONE;
}

//- map behavior as primitives ------------------------------------------------
PyObject* MapContains(PyObject* self, PyObject* obj)
{
// Implement python's __contains__ for std::map<>s
    PyObject* result = nullptr;

    PyObject* iter = CallPyObjMethod(self, "find", obj);
    if (CPPInstance_Check(iter)) {
        PyObject* end = PyObject_CallMethodObjArgs(self, PyStrings::gEnd, nullptr);
        if (CPPInstance_Check(end)) {
            if (!PyObject_RichCompareBool(iter, end, Py_EQ)) {
                Py_INCREF(Py_True);
                result = Py_True;
            }
        }
        Py_XDECREF(end);
    }
    Py_XDECREF(iter);

    if (!result) {
        PyErr_Clear();            // e.g. wrong argument type, which should always lead to False
        Py_INCREF(Py_False);
        result = Py_False;
    }

    return result;
}

//- STL container iterator support --------------------------------------------
PyObject* StlSequenceIter(PyObject* self)
{
// Implement python's __iter__ for std::iterator<>s
    PyObject* iter = PyObject_CallMethodObjArgs(self, PyStrings::gBegin, nullptr);
    if (iter) {
        PyObject* end = PyObject_CallMethodObjArgs(self, PyStrings::gEnd, nullptr);
        if (end)
            PyObject_SetAttr(iter, PyStrings::gEnd, end);
        Py_XDECREF(end);

    // add iterated collection as attribute so its refcount stays >= 1 while it's being iterated over
        PyObject_SetAttr(iter, CPyCppyy_PyText_FromString("_collection"), self);
    }
    return iter;
}

//- generic iterator support over a sequence with operator[] and size ---------
//-----------------------------------------------------------------------------
static PyObject* index_iter(PyObject* c) {
    indexiterobject* ii = PyObject_GC_New(indexiterobject, &IndexIter_Type);
    if (!ii) return nullptr;

    Py_INCREF(c);
    ii->ii_container = c;
    ii->ii_pos = 0;
    ii->ii_len = PySequence_Size(c);

    PyObject_GC_Track(ii);
    return (PyObject*)ii;
}


//- safe indexing for STL-like vector w/o iterator dictionaries ---------------
/* replaced by indexiterobject iteration, but may still have some future use ...
PyObject* CheckedGetItem(PyObject* self, PyObject* obj)
{
// Implement a generic python __getitem__ for STL-like classes that are missing the
// reflection info for their iterators. This is then used for iteration by means of
// consecutive indeces, it such index is of integer type.
    Py_ssize_t size = PySequence_Size(self);
    Py_ssize_t idx  = PyInt_AsSsize_t(obj);
    if ((size == (Py_ssize_t)-1 || idx == (Py_ssize_t)-1) && PyErr_Occurred()) {
    // argument conversion problem: let method itself resolve anew and report
        PyErr_Clear();
        return PyObject_CallMethodObjArgs(self, PyStrings::gGetNoCheck, obj, nullptr);
    }

    bool inbounds = false;
    if (idx < 0) idx += size;
    if (0 <= idx && 0 <= size && idx < size)
        inbounds = true;

    if (inbounds)
        return PyObject_CallMethodObjArgs(self, PyStrings::gGetNoCheck, obj, nullptr);
    else
        PyErr_SetString( PyExc_IndexError, "index out of range" );

    return nullptr;
}*/

//- pair as sequence to allow tuple unpacking --------------------------------
PyObject* PairUnpack(PyObject* self, PyObject* pyindex)
{
// For std::map<> iteration, unpack std::pair<>s into tuples for the loop.
    long idx = PyLong_AsLong(pyindex);
    if (idx == -1 && PyErr_Occurred())
        return nullptr;

    if (!CPPInstance_Check(self) || !((CPPInstance*)self)->GetObject()) {
        PyErr_SetString(PyExc_TypeError, "unsubscriptable object");
        return nullptr;
    }

    if ((int)idx == 0)
        return PyObject_GetAttr(self, PyStrings::gFirst);
    else if ((int)idx == 1)
        return PyObject_GetAttr(self, PyStrings::gSecond);

// still here? Trigger stop iteration
    PyErr_SetString(PyExc_IndexError, "out of bounds");
    return nullptr;
}

//- simplistic len() functions -----------------------------------------------
PyObject* ReturnTwo(CPPInstance*, PyObject*) {
    return PyInt_FromLong(2);
}


//- shared_ptr behavior --------------------------------------------------------
PyObject* SharedPtrInit(PyObject* self, PyObject* args, PyObject* /* kwds */)
{
// since the shared pointer will take ownership, we need to relinquish it
    PyObject* realInit = PyObject_GetAttrString(self, "__real_init");
    if (realInit) {
        PyObject* result = PyObject_Call(realInit, args, nullptr);
        Py_DECREF(realInit);
        if (result && PyTuple_GET_SIZE(args) == 1 && CPPInstance_Check(PyTuple_GET_ITEM(args, 0)))
            PyObject_SetAttrString(PyTuple_GET_ITEM(args, 0), "__python_owns__", Py_False);
        return result;
    }
    return nullptr;
}


//- string behavior as primitives --------------------------------------------
#if PY_VERSION_HEX >= 0x03000000
// TODO: this is wrong, b/c it doesn't order
static int PyObject_Compare(PyObject* one, PyObject* other) {
    return !PyObject_RichCompareBool(one, other, Py_EQ);
}
#endif
static inline PyObject* CPyCppyy_PyString_FromCppString(std::string* s) {
    return CPyCppyy_PyText_FromStringAndSize(s->c_str(), s->size());
}

static inline PyObject* CPyCppyy_PyString_FromCppString(std::wstring* s) {
    return PyUnicode_FromWideChar(s->c_str(), s->size());
}

#define CPPYY_IMPL_STRING_PYTHONIZATION(type, name)                          \
static PyObject* name##StringGetData(PyObject* self)                         \
{                                                                            \
    if (CPyCppyy::CPPInstance_Check(self)) {                                 \
        type* obj = ((type*)((CPPInstance*)self)->GetObject());              \
        if (obj) {                                                           \
            return CPyCppyy_PyString_FromCppString(obj);                     \
        } else {                                                             \
            return CPPInstance_Type.tp_str(self);                            \
        }                                                                    \
    }                                                                        \
    PyErr_Format(PyExc_TypeError, "object mismatch (%s expected)", #type);   \
    return nullptr;                                                          \
}                                                                            \
                                                                             \
PyObject* name##StringRepr(PyObject* self)                                   \
{                                                                            \
    PyObject* data = name##StringGetData(self);                              \
    if (data) {                                                              \
        PyObject* repr = PyObject_Repr(data);                                \
        Py_DECREF(data);                                                     \
        return repr;                                                         \
    }                                                                        \
    return nullptr;                                                          \
}                                                                            \
                                                                             \
PyObject* name##StringIsEqual(PyObject* self, PyObject* obj)                 \
{                                                                            \
    PyObject* data = name##StringGetData(self);                              \
    if (data) {                                                              \
        PyObject* result = PyObject_RichCompare(data, obj, Py_EQ);           \
        Py_DECREF(data);                                                     \
        return result;                                                       \
    }                                                                        \
    return nullptr;                                                          \
}                                                                            \
                                                                             \
PyObject* name##StringIsNotEqual(PyObject* self, PyObject* obj)              \
{                                                                            \
    PyObject* data = name##StringGetData(self);                              \
    if (data) {                                                              \
        PyObject* result = PyObject_RichCompare(data, obj, Py_NE);           \
        Py_DECREF(data);                                                     \
        return result;                                                       \
    }                                                                        \
    return nullptr;                                                          \
}

// Only define StlStringCompare:
#define CPPYY_IMPL_STRING_PYTHONIZATION_CMP(type, name)                      \
CPPYY_IMPL_STRING_PYTHONIZATION(type, name)                                  \
PyObject* name##StringCompare(PyObject* self, PyObject* obj)                 \
{                                                                            \
    PyObject* data = name##StringGetData(self);                              \
    int result = 0;                                                          \
    if (data) {                                                              \
        result = PyObject_Compare(data, obj);                                \
        Py_DECREF(data);                                                     \
    }                                                                        \
    if (PyErr_Occurred())                                                    \
        return nullptr;                                                      \
    return PyInt_FromLong(result);                                           \
}

CPPYY_IMPL_STRING_PYTHONIZATION_CMP(std::string, Stl)
CPPYY_IMPL_STRING_PYTHONIZATION_CMP(std::wstring, StlW)

Py_hash_t StlStringHash(PyObject* self)
{
// std::string objects hash to the same values as Python strings to allow
// matches in dictionaries etc.
    PyObject* data = StlStringGetData(self);
    Py_hash_t h = CPyCppyy_PyText_Type.tp_hash(data);
    Py_DECREF(data);
    return h;
}


//- STL iterator behavior ----------------------------------------------------
PyObject* StlIterNext(PyObject* self)
{
// Python iterator protocol __next__ for STL forward iterators.
    PyObject* next = nullptr;
    PyObject* last = PyObject_GetAttr(self, PyStrings::gEnd);

    if (last) {
    // handle special case of empty container (i.e. self is end)
        if (PyObject_RichCompareBool(last, self, Py_EQ) == 0) {
        // first, get next from the _current_ iterator as internal state may change
        // when call post or pre increment
            next = PyObject_CallMethodObjArgs(self, PyStrings::gDeref, nullptr);
            if (!next) PyErr_Clear();

        // use postinc, even as the C++11 range-based for loops prefer preinc b/c
        // that allows the current value from the iterator to be had from __deref__,
        // an issue that does not come up in C++
            static PyObject* dummy = PyInt_FromLong(1l);
            PyObject* iter = PyObject_CallMethodObjArgs(self, PyStrings::gPostInc, dummy, nullptr);
            if (!iter) {
            // allow preinc, as in that case likely __deref__ is not defined and it
            // is the iterator rather that is returned in the loop
                PyErr_Clear();
                iter = PyObject_CallMethodObjArgs(self, PyStrings::gPreInc, nullptr);
            }
            if (iter) {
            // prefer != as per C++11 range-based for
                int isNotEnd = PyObject_RichCompareBool(last, iter, Py_NE);
                if (isNotEnd && !next) {
                // if no dereference, continue iterating over the iterator
                    Py_INCREF(iter);
                    next = iter;
                }
                Py_DECREF(iter);
            } else {
            // fail current next, even if available
                Py_XDECREF(next);
                next = nullptr;
            }
        } else {
            PyErr_SetString(PyExc_StopIteration, "");
        }
        Py_DECREF(last);
    }

    if (!next) PyErr_SetString(PyExc_StopIteration, "");
    return next;
}


//- STL complex<T> behavior --------------------------------------------------
#define COMPLEX_METH_GETSET(name, cppname)                                   \
static PyObject* name##ComplexGet(PyObject* self, void*) {                   \
    return PyObject_CallMethodObjArgs(self, cppname, nullptr);               \
}                                                                            \
static int name##ComplexSet(PyObject* self, PyObject* value, void*) {        \
    PyObject* result = PyObject_CallMethodObjArgs(self, cppname, value, nullptr);\
    if (result) {                                                            \
        Py_DECREF(result);                                                   \
        return 0;                                                            \
    }                                                                        \
    return -1;                                                               \
}                                                                            \
PyGetSetDef name##Complex{(char*)#name, (getter)name##ComplexGet, (setter)name##ComplexSet, nullptr, nullptr};

COMPLEX_METH_GETSET(real, PyStrings::gCppReal)
COMPLEX_METH_GETSET(imag, PyStrings::gCppImag)

static PyObject* ComplexComplex(PyObject* self) {
    PyObject* real = PyObject_CallMethodObjArgs(self, PyStrings::gCppReal, nullptr);
    if (!real) return nullptr;
    double r = PyFloat_AsDouble(real);
    Py_DECREF(real);
    if (r == -1. && PyErr_Occurred())
        return nullptr;

    PyObject* imag = PyObject_CallMethodObjArgs(self, PyStrings::gCppImag, nullptr);
    if (!imag) return nullptr;
    double i = PyFloat_AsDouble(imag);
    Py_DECREF(imag);
    if (i == -1. && PyErr_Occurred())
        return nullptr;

    return PyComplex_FromDoubles(r, i);
}

static PyObject* ComplexRepr(PyObject* self) {
    PyObject* real = PyObject_CallMethodObjArgs(self, PyStrings::gCppReal, nullptr);
    if (!real) return nullptr;
    double r = PyFloat_AsDouble(real);
    Py_DECREF(real);
    if (r == -1. && PyErr_Occurred())
        return nullptr;

    PyObject* imag = PyObject_CallMethodObjArgs(self, PyStrings::gCppImag, nullptr);
    if (!imag) return nullptr;
    double i = PyFloat_AsDouble(imag);
    Py_DECREF(imag);
    if (i == -1. && PyErr_Occurred())
        return nullptr;

    std::ostringstream s;
    s << '(' << r << '+' << i << "j)";
    return CPyCppyy_PyText_FromString(s.str().c_str());
}

static PyObject* ComplexDRealGet(CPPInstance* self, void*)
{
    return PyFloat_FromDouble(((std::complex<double>*)self->GetObject())->real());
}

static int ComplexDRealSet(CPPInstance* self, PyObject* value, void*)
{
    double d = PyFloat_AsDouble(value);
    if (d == -1.0 && PyErr_Occurred())
        return -1;
    ((std::complex<double>*)self->GetObject())->real(d);
    return 0;
}

PyGetSetDef ComplexDReal{(char*)"real", (getter)ComplexDRealGet, (setter)ComplexDRealSet, nullptr, nullptr};


static PyObject* ComplexDImagGet(CPPInstance* self, void*)
{
    return PyFloat_FromDouble(((std::complex<double>*)self->GetObject())->imag());
}

static int ComplexDImagSet(CPPInstance* self, PyObject* value, void*)
{
    double d = PyFloat_AsDouble(value);
    if (d == -1.0 && PyErr_Occurred())
        return -1;
    ((std::complex<double>*)self->GetObject())->imag(d);
    return 0;
}

PyGetSetDef ComplexDImag{(char*)"imag", (getter)ComplexDImagGet, (setter)ComplexDImagSet, nullptr, nullptr};

static PyObject* ComplexDComplex(CPPInstance* self)
{
    double r = ((std::complex<double>*)self->GetObject())->real();
    double i = ((std::complex<double>*)self->GetObject())->imag();
    return PyComplex_FromDoubles(r, i);
}


} // unnamed namespace


//- public functions ---------------------------------------------------------
namespace CPyCppyy {
    std::set<std::string> gIteratorTypes;
}

bool CPyCppyy::Pythonize(PyObject* pyclass, const std::string& name)
{
// Add pre-defined pythonizations (for STL and ROOT) to classes based on their
// signature and/or class name.
    if (!pyclass)
        return false;

   CPPScope* klass = (CPPScope*)pyclass;

//- method name based pythonization ------------------------------------------

// for smart pointer style classes that are otherwise not known as such; would
// prefer operator-> as that returns a pointer (which is simpler since it never
// has to deal with ref-assignment), but operator* plays better with STL iters
// and algorithms
    if (HasAttrDirect(pyclass, PyStrings::gDeref) && !Cppyy::IsSmartPtr(klass->fCppType))
        Utility::AddToClass(pyclass, "__getattr__", (PyCFunction)DeRefGetAttr, METH_O);
    else if (HasAttrDirect(pyclass, PyStrings::gFollow) && !Cppyy::IsSmartPtr(klass->fCppType))
        Utility::AddToClass(pyclass, "__getattr__", (PyCFunction)FollowGetAttr, METH_O);

// for STL containers, and user classes modeled after them
    if (HasAttrDirect(pyclass, PyStrings::gSize))
        Utility::AddToClass(pyclass, "__len__", "size");

    if (!IsTemplatedSTLClass(name, "vector")  &&      // vector is dealt with below
           !((PyTypeObject*)pyclass)->tp_iter) {
        if (HasAttrDirect(pyclass, PyStrings::gBegin) && HasAttrDirect(pyclass, PyStrings::gEnd)) {
        // obtain the name of the return type
            const auto& v = Cppyy::GetMethodIndicesFromName(klass->fCppType, "begin");
            if (!v.empty()) {
            // check return type; if not explicitly an iterator, add it to the "known" return
            // types to add the "next" method on use
                Cppyy::TCppMethod_t meth = Cppyy::GetMethod(klass->fCppType, v[0]);
                const std::string& resname = Cppyy::GetMethodResultType(meth);
                if (Cppyy::GetScope(resname)) {
                    if (resname.find("iterator") == std::string::npos)
                        gIteratorTypes.insert(resname);

                // install iterator protocol a la STL
                    ((PyTypeObject*)pyclass)->tp_iter = (getiterfunc)StlSequenceIter;
                    Utility::AddToClass(pyclass, "__iter__", (PyCFunction)StlSequenceIter, METH_NOARGS);
                }
            }
        }
        if (!((PyTypeObject*)pyclass)->tp_iter &&     // no iterator resolved
                HasAttrDirect(pyclass, PyStrings::gGetItem) && HasAttrDirect(pyclass, PyStrings::gLen)) {
        // Python will iterate over __getitem__ using integers, but C++ operator[] will never raise
        // a StopIteration. A checked getitem (raising IndexError if beyond size()) works in some
        // cases but would mess up if operator[] is meant to implement an associative container. So,
        // this has to be implemented as an interator protocol.
            ((PyTypeObject*)pyclass)->tp_iter = (getiterfunc)index_iter;
            Utility::AddToClass(pyclass, "__iter__", (PyCFunction)index_iter, METH_NOARGS);
        }
    }

// operator==/!= are used in op_richcompare of CPPInstance, which subsequently allows
// comparisons to None; if no operator is available, a hook is installed for lazy
// lookups in the global and/or class namespace
    if (HasAttrDirect(pyclass, PyStrings::gEq, true)) {
        PyObject* cppol = PyObject_GetAttr(pyclass, PyStrings::gEq);
        if (!klass->fOperators) klass->fOperators = new Utility::PyOperators();
        klass->fOperators->fEq = cppol;
    // re-insert the forwarding __eq__ from the CPPInstance in case there was a Python-side
    // override in the base class
        static PyObject* top_eq = nullptr;
        if (!top_eq) {
            PyObject* top_cls = PyObject_GetAttrString(gThisModule, "CPPInstance");
            top_eq = PyObject_GetAttr(top_cls, PyStrings::gEq);
            Py_DECREF(top_eq);    // make it borrowed
            Py_DECREF(top_cls);
        }
        PyObject_SetAttr(pyclass, PyStrings::gEq, top_eq);
    }

    if (HasAttrDirect(pyclass, PyStrings::gNe, true)) {
        PyObject* cppol = PyObject_GetAttr(pyclass, PyStrings::gNe);
        if (!klass->fOperators) klass->fOperators = new Utility::PyOperators();
        klass->fOperators->fNe = cppol;
    // re-insert the forwarding __ne__ (same reason as above for __eq__)
        static PyObject* top_ne = nullptr;
        if (!top_ne) {
            PyObject* top_cls = PyObject_GetAttrString(gThisModule, "CPPInstance");
            top_ne = PyObject_GetAttr(top_cls, PyStrings::gNe);
            Py_DECREF(top_ne);    // make it borrowed
            Py_DECREF(top_cls);
        }
        PyObject_SetAttr(pyclass, PyStrings::gNe, top_ne);
    }


//- class name based pythonization -------------------------------------------

    if (IsTemplatedSTLClass(name, "vector")) {

    // std::vector<bool> is a special case in C++
        if (!sVectorBoolTypeID) sVectorBoolTypeID = (Cppyy::TCppType_t)Cppyy::GetScope("std::vector<bool>");
        if (klass->fCppType == sVectorBoolTypeID) {
            Utility::AddToClass(pyclass, "__getitem__", (PyCFunction)VectorBoolGetItem, METH_O);
            Utility::AddToClass(pyclass, "__setitem__", (PyCFunction)VectorBoolSetItem);
        } else {
        // constructor that takes python collections
            Utility::AddToClass(pyclass, "__real_init", "__init__");
            Utility::AddToClass(pyclass, "__init__", (PyCFunction)VectorInit, METH_VARARGS | METH_KEYWORDS);

        // data with size
            Utility::AddToClass(pyclass, "__real_data", "data");
            Utility::AddToClass(pyclass, "data", (PyCFunction)VectorData);

        // checked getitem
            if (HasAttrDirect(pyclass, PyStrings::gLen)) {
                Utility::AddToClass(pyclass, "_getitem__unchecked", "__getitem__");
                Utility::AddToClass(pyclass, "__getitem__", (PyCFunction)VectorGetItem, METH_O);
            }

        // vector-optimized iterator protocol
            ((PyTypeObject*)pyclass)->tp_iter      = (getiterfunc)vector_iter;

        // helpers for iteration
            const std::string& vtype = Cppyy::ResolveName(name+"::value_type");
            size_t typesz = Cppyy::SizeOf(vtype);
            if (typesz) {
                PyObject* pyvalue_size = PyLong_FromSsize_t(typesz);
                PyObject_SetAttrString(pyclass, "value_size", pyvalue_size);
                Py_DECREF(pyvalue_size);

                PyObject* pyvalue_type = CPyCppyy_PyText_FromString(vtype.c_str());
                PyObject_SetAttrString(pyclass, "value_type", pyvalue_type);
                Py_DECREF(pyvalue_type);
            }
        }
    }

    else if (IsTemplatedSTLClass(name, "map")) {
        Utility::AddToClass(pyclass, "__contains__", (PyCFunction)MapContains, METH_O);
    }

    else if (IsTemplatedSTLClass(name, "pair")) {
        Utility::AddToClass(pyclass, "__getitem__", (PyCFunction)PairUnpack, METH_O);
        Utility::AddToClass(pyclass, "__len__", (PyCFunction)ReturnTwo, METH_NOARGS);
    }

    if (IsTemplatedSTLClass(name, "shared_ptr")) {
        Utility::AddToClass(pyclass, "__real_init", "__init__");
        Utility::AddToClass(pyclass, "__init__", (PyCFunction)SharedPtrInit, METH_VARARGS | METH_KEYWORDS);
    }

    else if (name.find("iterator") != std::string::npos || gIteratorTypes.find(name) != gIteratorTypes.end()) {
        ((PyTypeObject*)pyclass)->tp_iternext = (iternextfunc)StlIterNext;
        Utility::AddToClass(pyclass, CPPYY__next__, (PyCFunction)StlIterNext, METH_NOARGS);
        ((PyTypeObject*)pyclass)->tp_iter = (getiterfunc)PyObject_SelfIter;
        Utility::AddToClass(pyclass, "__iter__", (PyCFunction)PyObject_SelfIter, METH_NOARGS);
    }

    else if (name == "string" || name == "std::string") { // TODO: ask backend as well
        Utility::AddToClass(pyclass, "__repr__", (PyCFunction)StlStringRepr,    METH_NOARGS);
        Utility::AddToClass(pyclass, "__str__",  (PyCFunction)StlStringGetData, METH_NOARGS);
        Utility::AddToClass(pyclass, "__cmp__",  (PyCFunction)StlStringCompare, METH_O);
        Utility::AddToClass(pyclass, "__eq__",   (PyCFunction)StlStringIsEqual, METH_O);
        Utility::AddToClass(pyclass, "__ne__",   (PyCFunction)StlStringIsNotEqual, METH_O);
        ((PyTypeObject*)pyclass)->tp_hash = (hashfunc)StlStringHash;
    }

    else if (name == "basic_string<wchar_t,char_traits<wchar_t>,allocator<wchar_t> >" || \
             name == "std::basic_string<wchar_t,char_traits<wchar_t>,allocator<wchar_t> >") {
        Utility::AddToClass(pyclass, "__repr__", (PyCFunction)StlWStringRepr,    METH_NOARGS);
        Utility::AddToClass(pyclass, "__str__",  (PyCFunction)StlWStringGetData, METH_NOARGS);
        Utility::AddToClass(pyclass, "__cmp__",  (PyCFunction)StlWStringCompare, METH_O);
        Utility::AddToClass(pyclass, "__eq__",   (PyCFunction)StlWStringIsEqual, METH_O);
        Utility::AddToClass(pyclass, "__ne__",   (PyCFunction)StlWStringIsNotEqual, METH_O);
    }

    else if (name == "complex<double>" || name == "std::complex<double>") {
        Utility::AddToClass(pyclass, "__cpp_real", "real");
        PyObject_SetAttrString(pyclass, "real",  PyDescr_NewGetSet((PyTypeObject*)pyclass, &ComplexDReal));
        Utility::AddToClass(pyclass, "__cpp_imag", "imag");
        PyObject_SetAttrString(pyclass, "imag",  PyDescr_NewGetSet((PyTypeObject*)pyclass, &ComplexDImag));
        Utility::AddToClass(pyclass, "__complex__", (PyCFunction)ComplexDComplex, METH_NOARGS);
        Utility::AddToClass(pyclass, "__repr__", (PyCFunction)ComplexRepr, METH_NOARGS);
    }

    else if (IsTemplatedSTLClass(name, "complex")) {
        Utility::AddToClass(pyclass, "__cpp_real", "real");
        PyObject_SetAttrString(pyclass, "real", PyDescr_NewGetSet((PyTypeObject*)pyclass, &realComplex));
        Utility::AddToClass(pyclass, "__cpp_imag", "imag");
        PyObject_SetAttrString(pyclass, "imag", PyDescr_NewGetSet((PyTypeObject*)pyclass, &imagComplex)); 
        Utility::AddToClass(pyclass, "__complex__", (PyCFunction)ComplexComplex, METH_NOARGS);
        Utility::AddToClass(pyclass, "__repr__", (PyCFunction)ComplexRepr, METH_NOARGS);
    }

// direct user access; there are two calls here:
//   - explicit pythonization: won't fall through to the base classes and is preferred if present
//   - normal pythonization: only called if explicit isn't present, falls through to base classes
    bool bUserOk = true; PyObject* res = nullptr;
    PyObject* pyname = CPyCppyy_PyText_FromString(name.c_str());
    if (HasAttrDirect(pyclass, PyStrings::gExPythonize)) {
        res = PyObject_CallMethodObjArgs(pyclass, PyStrings::gExPythonize, pyclass, pyname, nullptr);
        bUserOk = (bool)res;
    } else {
        PyObject* func = PyObject_GetAttr(pyclass, PyStrings::gPythonize);
        if (func) {
            res = PyObject_CallFunctionObjArgs(func, pyclass, pyname, nullptr);
            Py_DECREF(func);
            bUserOk = (bool)res;
        } else
            PyErr_Clear();
    }
    if (!bUserOk) {
        Py_DECREF(pyname);
        return false;
    } else {
        Py_XDECREF(res);
        // pyname handed to tuple below
    }

// call registered pythonizors, if any
    PyObject* args = PyTuple_New(2);
    Py_INCREF(pyclass);
    PyTuple_SET_ITEM(args, 0, pyclass);

    std::string outer_scope = TypeManip::extract_namespace(name);

    bool pstatus = true;
    auto p = outer_scope.empty() ? gPythonizations.end() : gPythonizations.find(outer_scope);
    if (p == gPythonizations.end()) {
        p = gPythonizations.find("");
        PyTuple_SET_ITEM(args, 1, pyname);
    } else {
        PyTuple_SET_ITEM(args, 1, CPyCppyy_PyText_FromString(
                                      name.substr(outer_scope.size()+2, std::string::npos).c_str()));
        Py_DECREF(pyname);
    }

    if (p != gPythonizations.end()) {
        for (auto pythonizor : p->second) {
            PyObject* result = PyObject_CallObject(pythonizor, args);
            if (!result) {
            // TODO: detail error handling for the pythonizors
                pstatus = false;
                break;
            }
            Py_DECREF(result);
        }
    }

    Py_DECREF(args);

// phew! all done ...
    return pstatus;
}
