// Bindings
#include "CPyCppyy.h"
#include "LowLevelViews.h"
#include "Converters.h"

// Standard
#include <map>
#include <assert.h>
#include <limits.h>


//= memoryview-like object ===================================================
// This is largely setup just like Python builtin memory view objects, with
// the exceptions that there is no need of a "base" object (it views on C++
// memory, not a Python object with a buffer interface), it uses the CPyCppyy
// converters, and typed results and assignments are supported. Code reused
// under PSF License Version 2.

//= CPyCppyy low level view construction/destruction =========================
static CPyCppyy::LowLevelView* ll_new(PyTypeObject* subtype, PyObject*, PyObject*)
{
// Create a new low level ptr type
    CPyCppyy::LowLevelView* pyobj = (CPyCppyy::LowLevelView*)subtype->tp_alloc(subtype, 0);
    if (!pyobj) PyErr_Print();
    memset(&pyobj->fBufInfo, 0, sizeof(Py_buffer));
    pyobj->fBuf = nullptr;
    pyobj->fConverter = nullptr;

    return pyobj;
}

//----------------------------------------------------------------------------
static void ll_dealloc(CPyCppyy::LowLevelView* pyobj)
{
// Destruction requires the deletion of the converter (if any)
    PyMem_Free(pyobj->fBufInfo.shape);
    PyMem_Free(pyobj->fBufInfo.strides);
    if (pyobj->fConverter && pyobj->fConverter->HasState()) delete pyobj->fConverter;
    Py_TYPE(pyobj)->tp_free((PyObject*)pyobj);
}


//---------------------------------------------------------------------------
static PyObject* ll_typecode(CPyCppyy::LowLevelView* self, void*)
{
    return CPyCppyy_PyText_FromString((char*)self->fBufInfo.format);
}

//---------------------------------------------------------------------------
static PyGetSetDef ll_getset[] = {
    {(char*)"format",   (getter)ll_typecode, nullptr, nullptr, nullptr},
    {(char*)"typecode", (getter)ll_typecode, nullptr, nullptr, nullptr},
    {(char*)nullptr, nullptr, nullptr, nullptr, nullptr }
};


//---------------------------------------------------------------------------
static PyObject* ll_reshape(CPyCppyy::LowLevelView* self, PyObject* shape)
{
// Allow the user to fix up the actual (type-strided) size of the buffer.
    if (!PyTuple_Check(shape) || PyTuple_GET_SIZE(shape) != 1) {
        if (shape) {
            PyObject* pystr = PyObject_Str(shape);
            if (pystr) {
                PyErr_Format(PyExc_TypeError, "tuple object of length 1 expected, received %s",
                    CPyCppyy_PyText_AsStringChecked(pystr));
                Py_DECREF(pystr);
                return nullptr;
            }
        }
        PyErr_SetString(PyExc_TypeError, "tuple object of length 1 expected");
        return nullptr;
    }

    Py_ssize_t nlen = PyInt_AsSsize_t(PyTuple_GET_ITEM(shape, 0));
    if (nlen == -1 && PyErr_Occurred())
        return nullptr;

    self->fBufInfo.len = nlen * self->fBufInfo.itemsize;
    if (self->fBufInfo.ndim == 1 && self->fBufInfo.shape)
        self->fBufInfo.shape[0] = nlen;
    else {
        PyErr_SetString(PyExc_TypeError, "unsupported buffer dimensions");
        return nullptr;
    }

    Py_RETURN_NONE;
}

//---------------------------------------------------------------------------
static PyMethodDef ll_methods[] = {
    {(char*)"reshape",     (PyCFunction)ll_reshape, METH_O,      nullptr},
    {(char*)nullptr, nullptr, 0, nullptr}
};


//- Copy memoryview buffers =================================================

// The functions in this section take a source and a destination buffer
// with the same logical structure: format, itemsize, ndim and shape
// are identical, with ndim > 0.

// Check for the presence of suboffsets in the first dimension.
#define HAVE_PTR(suboffsets, dim) (suboffsets && suboffsets[dim] >= 0)
// Adjust ptr if suboffsets are present.
#define ADJUST_PTR(ptr, suboffsets, dim) \
    (HAVE_PTR(suboffsets, dim) ? *((char**)ptr) + suboffsets[dim] : ptr)

// Assumptions: ndim >= 1. The macro tests for a corner case that should
// perhaps be explicitly forbidden in the PEP.
#define HAVE_SUBOFFSETS_IN_LAST_DIM(view) \
    (view->suboffsets && view->suboffsets[dest->ndim-1] >= 0)

//---------------------------------------------------------------------------
static inline int last_dim_is_contiguous(const Py_buffer *dest, const Py_buffer *src)
{
    assert(dest->ndim > 0 && src->ndim > 0);
    return (!HAVE_SUBOFFSETS_IN_LAST_DIM(dest) &&
            !HAVE_SUBOFFSETS_IN_LAST_DIM(src) &&
            dest->strides[dest->ndim-1] == dest->itemsize &&
            src->strides[src->ndim-1] == src->itemsize);
}

//---------------------------------------------------------------------------
static inline bool equiv_shape(const Py_buffer* dest, const Py_buffer* src)
{
// Two shapes are equivalent if they are either equal or identical up
// to a zero element at the same position. For example, in NumPy arrays
// the shapes [1, 0, 5] and [1, 0, 7] are equivalent.
    if (dest->ndim != src->ndim)
        return false;

    for (int i = 0; i < dest->ndim; i++) {
        if (dest->shape[i] != src->shape[i])
            return 0;
        if (dest->shape[i] == 0)
            break;
    }

    return true;
}

//---------------------------------------------------------------------------
static bool equiv_structure(const Py_buffer* dest, const Py_buffer* src)
{
// Check that the logical structure of the destination and source buffers
// is identical.
    if (strcmp(dest->format, src->format) != 0 || dest->itemsize != src->itemsize ||
        !equiv_shape(dest, src)) {
        PyErr_SetString(PyExc_ValueError,
            "low level pointer assignment: lvalue and rvalue have different structures");
        return false;
    }

    return true;
}

//---------------------------------------------------------------------------
static void copy_base(const Py_ssize_t* shape, Py_ssize_t itemsize,
    char* dptr, const Py_ssize_t* dstrides, const Py_ssize_t* dsuboffsets,
    char* sptr, const Py_ssize_t* sstrides, const Py_ssize_t* ssuboffsets,
    char* mem)
{
// Base case for recursive multi-dimensional copying. Contiguous arrays are
// copied with very little overhead. Assumptions: ndim == 1, mem == nullptr or
// sizeof(mem) == shape[0] * itemsize.
    if (!mem) { // contiguous
        Py_ssize_t size = shape[0] * itemsize;
        if (dptr + size < sptr || sptr + size < dptr)
            memcpy(dptr, sptr, size); // no overlapping
        else
            memmove(dptr, sptr, size);
    }
    else {
        char *p;
        Py_ssize_t i;
        for (i=0, p=mem; i < shape[0]; p+=itemsize, sptr+=sstrides[0], i++) {
            char *xsptr = ADJUST_PTR(sptr, ssuboffsets, 0);
            memcpy(p, xsptr, itemsize);
        }
        for (i=0, p=mem; i < shape[0]; p+=itemsize, dptr+=dstrides[0], i++) {
            char *xdptr = ADJUST_PTR(dptr, dsuboffsets, 0);
            memcpy(xdptr, p, itemsize);
        }
    }

}

//---------------------------------------------------------------------------
static int copy_single(Py_buffer* dest, Py_buffer* src)
{
// Faster copying of one-dimensional arrays.
    char* mem = nullptr;

    assert(dest->ndim == 1);

    if (!equiv_structure(dest, src))
        return -1;

    if (!last_dim_is_contiguous(dest, src)) {
        mem = (char*)PyMem_Malloc(dest->shape[0] * dest->itemsize);
        if (!mem) {
            PyErr_NoMemory();
            return -1;
        }
    }

    copy_base(dest->shape, dest->itemsize,
              (char*)dest->buf, dest->strides, dest->suboffsets,
              (char*)src->buf, src->strides, src->suboffsets,
              mem);

    if (mem)
        PyMem_Free(mem);

    return 0;
}


//- Indexing and slicing ----------------------------------------------------
static char* lookup_dimension(Py_buffer& view, char* ptr, int dim, Py_ssize_t index)
{
    Py_ssize_t nitems; // items in the given dimension

    assert(view.shape);
    assert(view.strides);

    nitems = view.shape[dim];
    if (index < 0)
        index += nitems;

    if (index < 0 || index >= nitems) {
        PyErr_Format(PyExc_IndexError,
            "index out of bounds on dimension %d", dim + 1);
        return nullptr;
    }

    ptr += view.strides[dim] * index;
    ptr = ADJUST_PTR(ptr, view.suboffsets, dim);

    return ptr;
}

// Get the pointer to the item at index.
//---------------------------------------------------------------------------
static inline void* ptr_from_index(CPyCppyy::LowLevelView* llview, Py_ssize_t index)
{
    Py_buffer& view = llview->fBufInfo;
    return lookup_dimension(view, (char*)llview->get_buf(), 0, index);
}

// Get the pointer to the item at tuple.
//---------------------------------------------------------------------------
static void* ptr_from_tuple(CPyCppyy::LowLevelView* llview, PyObject* tup)
{
    Py_buffer& view = llview->fBufInfo;
    Py_ssize_t nindices = PyTuple_GET_SIZE(tup);
    if (nindices > view.ndim) {
        PyErr_Format(PyExc_TypeError,
            "cannot index %d-dimension view with %zd-element tuple", view.ndim, nindices);
        return nullptr;
    }

    char* ptr = (char*)llview->get_buf();
    for (Py_ssize_t dim = 0; dim < nindices; dim++) {
        Py_ssize_t index;
        index = PyNumber_AsSsize_t(PyTuple_GET_ITEM(tup, dim),
                                   PyExc_IndexError);
        if (index == -1 && PyErr_Occurred())
            return nullptr;
        ptr = lookup_dimension(view, ptr, (int)dim, index);
        if (!ptr)
            return nullptr;
    }
    return ptr;
}


//= mapping methods =========================================================
static Py_ssize_t ll_length(CPyCppyy::LowLevelView* self)
{
    if (!self->get_buf())
        return 0;
    return self->fBufInfo.ndim == 0 ? 1 : self->fBufInfo.shape[0];
}

//---------------------------------------------------------------------------
static inline int init_slice(Py_buffer* base, PyObject* _key, int dim)
{
    Py_ssize_t start, stop, step, slicelength;

#if PY_VERSION_HEX < 0x03000000
    PySliceObject* key = (PySliceObject*)_key;
#else
    PyObject* key = _key;
#endif

    if (PySlice_GetIndicesEx(key, base->shape[dim], &start, &stop, &step, &slicelength) < 0)
        return -1;

    if (!base->suboffsets || dim == 0) {
    adjust_buf:
        base->buf = (char *)base->buf + base->strides[dim] * start;
    }
    else {
        Py_ssize_t n = dim-1;
        while (n >= 0 && base->suboffsets[n] < 0)
            n--;
        if (n < 0)
            goto adjust_buf; // all suboffsets are negative
        base->suboffsets[n] = base->suboffsets[n] + base->strides[dim] * start;
    }
    base->shape[dim] = slicelength;
    base->strides[dim] = base->strides[dim] * step;

    return 0;
}

//---------------------------------------------------------------------------
static bool is_multislice(PyObject* key)
{
    if (!PyTuple_Check(key))
        return false;

    Py_ssize_t size = PyTuple_GET_SIZE(key);
    if (size == 0)
        return false;

    for (Py_ssize_t i = 0; i < size; i++) {
        PyObject *x = PyTuple_GET_ITEM(key, i);
        if (!PySlice_Check(x))
            return false;
    }
    return true;
}

//---------------------------------------------------------------------------
static Py_ssize_t is_multiindex(PyObject* key)
{
    if (!PyTuple_Check(key))
        return 0;

    Py_ssize_t size = PyTuple_GET_SIZE(key);
    for (Py_ssize_t i = 0; i < size; i++) {
        PyObject *x = PyTuple_GET_ITEM(key, i);
        if (!PyIndex_Check(x))
            return 0;
    }
    return 1;
}


// Return the item at index. In a one-dimensional view, this is an object
// with the type specified by view->format. Otherwise, the item is a sub-view.
// The function is used in ll_subscript() and ll_as_sequence.
//---------------------------------------------------------------------------
static PyObject* ll_item(CPyCppyy::LowLevelView* self, Py_ssize_t index)
{
    Py_buffer& view = self->fBufInfo;

    if (!self->get_buf()) {
        PyErr_SetString(PyExc_ReferenceError, "attempt to access a null-pointer");
        return nullptr;
    }

    if (view.ndim == 0) {
        PyErr_SetString(PyExc_TypeError, "invalid indexing of 0-dim memory");
        return nullptr;
    }

    void* ptr = ptr_from_index(self, index);
    if (ptr)
        return self->fConverter->FromMemory(ptr);

    return nullptr;      // error already set by lookup_dimension
}

// Return the item at position *key* (a tuple of indices).
//---------------------------------------------------------------------------
static PyObject* ll_item_multi(CPyCppyy::LowLevelView* self, PyObject *tup)
{
    Py_buffer& view = self->fBufInfo;
    Py_ssize_t nindices = PyTuple_GET_SIZE(tup);

    if (nindices < view.ndim) {
    // TODO: implement
        PyErr_SetString(PyExc_NotImplementedError,
                        "sub-views are not implemented");
        return nullptr;
    }

    void* ptr = ptr_from_tuple(self, tup);
    if (!ptr)
        return nullptr;
    return self->fConverter->FromMemory(ptr);
}


// llp[obj] returns an object holding the data for one element if obj
// fully indexes the lowlevelptr or another lowlevelptr object if it
// does not.
//
// 0-d lowlevelptr objects can be referenced using llp[...] or llp[()]
// but not with anything else.
//---------------------------------------------------------------------------
static PyObject* ll_subscript(CPyCppyy::LowLevelView* self, PyObject* key)
{
    Py_buffer& view = self->fBufInfo;

    if (view.ndim == 0) {
        if (PyTuple_Check(key) && PyTuple_GET_SIZE(key) == 0) {
            return self->fConverter->FromMemory(self->get_buf());
        }
        else if (key == Py_Ellipsis) {
            Py_INCREF(self);
            return (PyObject*)self;
        }
        else {
            PyErr_SetString(PyExc_TypeError,
                "invalid indexing of 0-dim memory");
            return nullptr;
        }
    }

    if (PyIndex_Check(key)) {
        Py_ssize_t index = PyNumber_AsSsize_t(key, PyExc_IndexError);
        if (index == -1 && PyErr_Occurred())
            return nullptr;
        return ll_item(self, index);
    }
    else if (PySlice_Check(key)) {
    // TODO: handle slicing. This should be simpler than the memoryview
    // case as there is no Python object holding the buffer.
        PyErr_SetString(PyExc_NotImplementedError,
            "multi-dimensional slicing is not implemented");
        return nullptr;
    }
    else if (is_multiindex(key)) {
        return ll_item_multi(self, key);
    }
    else if (is_multislice(key)) {
        PyErr_SetString(PyExc_NotImplementedError,
            "multi-dimensional slicing is not implemented");
        return nullptr;
    }

    PyErr_SetString(PyExc_TypeError, "invalid slice key");
    return nullptr;
}

//---------------------------------------------------------------------------
static int ll_ass_sub(CPyCppyy::LowLevelView* self, PyObject* key, PyObject* value)
{
    Py_buffer& view = self->fBufInfo;
    Py_buffer src;

    if (view.readonly) {
        PyErr_SetString(PyExc_TypeError, "cannot modify read-only memory");
        return -1;
    }

    if (value == nullptr) {
        PyErr_SetString(PyExc_TypeError, "cannot delete memory");
        return -1;
    }

    if (view.ndim == 0) {
        if (key == Py_Ellipsis ||
            (PyTuple_Check(key) && PyTuple_GET_SIZE(key) == 0)) {
            return self->fConverter->ToMemory(value, self->get_buf()) ? 0 : -1;
        }
        else {
            PyErr_SetString(PyExc_TypeError,
                "invalid indexing of 0-dim memory");
            return -1;
        }
    }

    if (PyIndex_Check(key)) {
        Py_ssize_t index;
        if (1 < view.ndim) {
            PyErr_SetString(PyExc_NotImplementedError,
                            "sub-views are not implemented");
            return -1;
        }
        index = PyNumber_AsSsize_t(key, PyExc_IndexError);
        if (index == -1 && PyErr_Occurred())
            return -1;
        void* ptr = ptr_from_index(self, index);
        if (ptr == nullptr)
            return -1;
        return self->fConverter->ToMemory(value, ptr) ? 0 : -1;
    }

    // one-dimensional: fast path
    if (PySlice_Check(key) && view.ndim == 1) {
        Py_buffer dest; // sliced view
        Py_ssize_t arrays[3];
        int ret = -1;

        // rvalue must be an exporter
        if (PyObject_GetBuffer(value, &src, PyBUF_FULL_RO) < 0)
            return ret;

        dest = view;
        dest.shape = &arrays[0]; dest.shape[0] = view.shape[0];
        dest.strides = &arrays[1]; dest.strides[0] = view.strides[0];
        if (view.suboffsets) {
            dest.suboffsets = &arrays[2]; dest.suboffsets[0] = view.suboffsets[0];
        }

        if (init_slice(&dest, key, 0) < 0)
            return -1;
        dest.len = dest.shape[0] * dest.itemsize;

        return copy_single(&dest, &src);
    }

    if (is_multiindex(key)) {
    // TODO: implement
        if (PyTuple_GET_SIZE(key) < view.ndim) {
            PyErr_SetString(PyExc_NotImplementedError,
                            "sub-views are not implemented");
            return -1;
        }
        void* ptr = ptr_from_tuple(self, key);
        if (ptr == nullptr)
            return -1;
        return self->fConverter->ToMemory(value, ptr) ? 0 : -1;
    }

    if (PySlice_Check(key) || is_multislice(key)) {
        // TODO: implement
        PyErr_SetString(PyExc_NotImplementedError,
            "LowLevelView slice assignments are currently restricted "
            "to ndim = 1");
        return -1;
    }

    PyErr_SetString(PyExc_TypeError, "invalid slice key");
    return -1;
}

#if PY_VERSION_HEX < 0x03000000
//---------------------------------------------------------------------------
static Py_ssize_t ll_oldgetbuf(CPyCppyy::LowLevelView* self, Py_ssize_t seg, void** pptr)
{
    if (seg != 0) {
        PyErr_SetString(PyExc_TypeError, "accessing non-existent segment");
        return -1;
    }

    *pptr = self->get_buf();
    return self->fBufInfo.len;
}

//---------------------------------------------------------------------------
static Py_ssize_t ll_getsegcount(PyObject*, Py_ssize_t* lenp)
{
    if (lenp) *lenp = 1;
    return 1;
}
#endif

//---------------------------------------------------------------------------
static int ll_getbuf(CPyCppyy::LowLevelView* self, Py_buffer* view, int flags)
{
// Simplified from memoryobject, as we're always dealing with C arrays.

// start with full copy
    *view = self->fBufInfo;

    if (!(flags & PyBUF_FORMAT)) {
        /* NULL indicates that the buffer's data type has been cast to 'B'.
           view->itemsize is the _previous_ itemsize. If shape is present,
           the equality product(shape) * itemsize = len still holds at this
           point. The equality calcsize(format) = itemsize does _not_ hold
           from here on! */
        view->format = NULL;
    }

    if ((flags & PyBUF_F_CONTIGUOUS) == PyBUF_F_CONTIGUOUS) {
        PyErr_SetString(PyExc_BufferError,
            "underlying buffer is not Fortran contiguous");
        return -1;
    }

    if (!(flags & PyBUF_FORMAT)) {
        /* PyBUF_SIMPLE or PyBUF_WRITABLE: at this point buf is C-contiguous,
           so base->buf = ndbuf->data. */
        if (view->format != NULL) {
            /* PyBUF_SIMPLE|PyBUF_FORMAT and PyBUF_WRITABLE|PyBUF_FORMAT do
               not make sense. */
            PyErr_Format(PyExc_BufferError,
                "cannot cast to unsigned bytes if the format flag is present");
            return -1;
        }
        /* product(shape) * itemsize = len and calcsize(format) = itemsize
           do _not_ hold from here on! */
        view->ndim = 1;
        view->shape = NULL;
    }

    view->obj = (PyObject*)self;
    Py_INCREF(view->obj);

    return 0;
}


//- mapping methods ---------------------------------------------------------
static PyMappingMethods ll_as_mapping = {
    (lenfunc)      ll_length,      // mp_length
    (binaryfunc)   ll_subscript,   // mp_subscript
    (objobjargproc)ll_ass_sub,     // mp_ass_subscript
};

//- sequence methods --------------------------------------------------------
static PySequenceMethods ll_as_sequence = {
    (lenfunc)ll_length,            // sq_length
    0,                             // sq_concat
    0,                             // sq_repeat
    (ssizeargfunc)ll_item,         // sq_item
    0,                             // sq_slice
    0,                             // sq_ass_item
    0,                             // sq_ass_slice
    0,                             // sq_contains
    0,                             // sq_inplace_concat
    0,                             // sq_inplace_repeat
};

//- buffer methods ----------------------------------------------------------
static PyBufferProcs ll_as_buffer = {
#if PY_VERSION_HEX < 0x03000000
    (readbufferproc)ll_oldgetbuf,   // bf_getreadbuffer
    (writebufferproc)ll_oldgetbuf,  // bf_getwritebuffer
    (segcountproc)ll_getsegcount,   // bf_getsegcount
    0,                              // bf_getcharbuffer
#endif
    (getbufferproc)ll_getbuf,       // bf_getbuffer
    0,                              // bf_releasebuffer
};


namespace CPyCppyy {

//= CPyCppyy low level view type ============================================
PyTypeObject LowLevelView_Type = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    (char*)"cppyy.LowLevelView",   // tp_name
    sizeof(CPyCppyy::LowLevelView),// tp_basicsize
    0,                             // tp_itemsize
    (destructor)ll_dealloc,        // tp_dealloc
    0,                             // tp_print
    0,                             // tp_getattr
    0,                             // tp_setattr
    0,                             // tp_compare
    0,                             // tp_repr
    0,                             // tp_as_number
    &ll_as_sequence,               // tp_as_sequence
    &ll_as_mapping,                // tp_as_mapping
    0,                             // tp_hash
    0,                             // tp_call
    0,                             // tp_str
    0,                             // tp_getattro
    0,                             // tp_setattro
    &ll_as_buffer,                 // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES |
        Py_TPFLAGS_BASETYPE,       // tp_flags
    (char*)"memory view on C++ pointer",     // tp_doc
    0,                             // tp_traverse
    0,                             // tp_clear
    0,                             // tp_richcompare
    0,                             // tp_weaklistoffset
    0,                             // tp_iter
    0,                             // tp_iternext
    ll_methods,                    // tp_methods
    0,                             // tp_members
    ll_getset,                     // tp_getset
    0,                             // tp_base
    0,                             // tp_dict
    0,                             // tp_descr_get
    0,                             // tp_descr_set
    0,                             // tp_dictoffset
    0,                             // tp_init
    0,                             // tp_alloc
    (newfunc)ll_new,               // tp_new
    0,                             // tp_free
    0,                             // tp_is_gc
    0,                             // tp_bases
    0,                             // tp_mro
    0,                             // tp_cache
    0,                             // tp_subclasses
    0                              // tp_weaklist
#if PY_VERSION_HEX >= 0x02030000
    , 0                            // tp_del
#endif
#if PY_VERSION_HEX >= 0x02060000
    , 0                            // tp_version_tag
#endif
#if PY_VERSION_HEX >= 0x03040000
    , 0                            // tp_finalize
#endif
};

} // namespace CPyCppyy

namespace {

template<typename T> struct typecode_traits {};
template<> struct typecode_traits<bool> {
    static constexpr const char* format = "?"; static constexpr const char* name = "bool"; };
template<> struct typecode_traits<signed char> {
    static constexpr const char* format = "b"; static constexpr const char* name = "signed char"; };
template<> struct typecode_traits<unsigned char> {
    static constexpr const char* format = "B"; static constexpr const char* name = "UCharAsInt"; };
#if __cplusplus > 201402L
template<> struct typecode_traits<std::byte> {
    static constexpr const char* format = "B"; static constexpr const char* name = "UCharAsInt"; };
#endif
template<> struct typecode_traits<const char*> {
    static constexpr const char* format = "b"; static constexpr const char* name = "const char*"; };
template<> struct typecode_traits<short> {
    static constexpr const char* format = "h"; static constexpr const char* name = "short"; };
template<> struct typecode_traits<unsigned short> {
    static constexpr const char* format = "H"; static constexpr const char* name = "unsigned short"; };
template<> struct typecode_traits<int> {
    static constexpr const char* format = "i"; static constexpr const char* name = "int"; };
template<> struct typecode_traits<unsigned int> {
    static constexpr const char* format = "I"; static constexpr const char* name = "unsigned int"; };
template<> struct typecode_traits<long> {
    static constexpr const char* format = "l"; static constexpr const char*
#if PY_VERSION_HEX < 0x03000000
        name = "int";
#else
        name = "long";
#endif
};
template<> struct typecode_traits<unsigned long> {
    static constexpr const char* format = "L"; static constexpr const char* name = "unsigned long"; };
template<> struct typecode_traits<long long> {
    static constexpr const char* format = "q"; static constexpr const char* name = "long long"; };
template<> struct typecode_traits<unsigned long long> {
    static constexpr const char* format = "Q"; static constexpr const char* name = "unsigned long long"; };
template<> struct typecode_traits<float> {
    static constexpr const char* format = "f"; static constexpr const char* name = "float"; };
template<> struct typecode_traits<double> {
    static constexpr const char* format = "d"; static constexpr const char* name = "double"; };
template<> struct typecode_traits<long double> {
    static constexpr const char* format = "D"; static constexpr const char* name = "long double"; };
template<> struct typecode_traits<std::complex<float>> {
    static constexpr const char* format = "Zf"; static constexpr const char* name = "std::complex<float>"; };
template<> struct typecode_traits<std::complex<double>> {
    static constexpr const char* format = "Zd"; static constexpr const char* name = "std::complex<double>"; };
template<> struct typecode_traits<std::complex<int>> {
    static constexpr const char* format = "Zi"; static constexpr const char* name = "std::complex<int>"; };
template<> struct typecode_traits<std::complex<long>> {
    static constexpr const char* format = "Zl"; static constexpr const char* name = "std::complex<long>"; };

} // unnamed namespace


//---------------------------------------------------------------------------
template<typename T>
static inline PyObject* CreateLowLevelViewT(T* address, Py_ssize_t* shape)
{
    using namespace CPyCppyy;
    Py_ssize_t nx = (shape && 0 <= shape[1]) ? shape[1] : INT_MAX/sizeof(T);
    PyObject* args = PyTuple_New(0);
    LowLevelView* llp =
        (LowLevelView*)LowLevelView_Type.tp_new(&LowLevelView_Type, args, nullptr);
    Py_DECREF(args);

    Py_buffer& view = llp->fBufInfo;
    view.buf            = address;
    view.obj            = nullptr;
    view.readonly       = 0;
    view.format         = (char*)typecode_traits<T>::format;
    view.ndim           = shape ? (int)shape[0] : 1;
    view.shape          = (Py_ssize_t*)PyMem_Malloc(view.ndim * sizeof(Py_ssize_t));
    view.shape[0]       = nx;      // view.len / view.itemsize
    view.strides        = (Py_ssize_t*)PyMem_Malloc(view.ndim * sizeof(Py_ssize_t));
    view.suboffsets     = nullptr;
    view.internal       = nullptr;

    if (view.ndim == 1) {
    // simple 1-dim array of the declared type
        view.len        = nx * sizeof(T);
        view.itemsize   = sizeof(T);
        llp->fConverter = CreateConverter(typecode_traits<T>::name);
    } else {
    // multi-dim array; sub-views are projected by using more LLViews
        view.len        = nx * sizeof(void*);
        view.itemsize   = sizeof(void*);

    // peel off one dimension and create a new LLView converter
        Py_ssize_t res = shape[1];
        shape[1] = shape[0] - 1;
        std::string tname{typecode_traits<T>::name};
        tname.append("*");        // make sure to ask for another array
    // TODO: although this will work, it means that "naive" loops are expensive
        llp->fConverter = CreateConverter(tname, &shape[1]);
        shape[1] = res;
    }

    view.strides[0]     = view.itemsize;

    return (PyObject*)llp;
}

//---------------------------------------------------------------------------
template<typename T>
static inline PyObject* CreateLowLevelViewT(T** address, Py_ssize_t* shape)
{
    using namespace CPyCppyy;
    T* buf = address ? *address : nullptr;
    LowLevelView* llp = (LowLevelView*)CreateLowLevelViewT(buf, shape);
    llp->set_buf((void**)address);
    return (PyObject*)llp;
}

//---------------------------------------------------------------------------
#define CPPYY_IMPL_VIEW_CREATOR(type)                                       \
PyObject* CPyCppyy::CreateLowLevelView(type* address, Py_ssize_t* shape) {  \
    return CreateLowLevelViewT<type>(address, shape);                       \
}                                                                           \
PyObject* CPyCppyy::CreateLowLevelView(type** address, Py_ssize_t* shape) { \
    return CreateLowLevelViewT<type>(address, shape);                       \
}

CPPYY_IMPL_VIEW_CREATOR(bool);
CPPYY_IMPL_VIEW_CREATOR(signed char);
CPPYY_IMPL_VIEW_CREATOR(unsigned char);
#if __cplusplus > 201402L
CPPYY_IMPL_VIEW_CREATOR(std::byte);
#endif
CPPYY_IMPL_VIEW_CREATOR(short);
CPPYY_IMPL_VIEW_CREATOR(unsigned short);
CPPYY_IMPL_VIEW_CREATOR(int);
CPPYY_IMPL_VIEW_CREATOR(unsigned int);
CPPYY_IMPL_VIEW_CREATOR(long);
CPPYY_IMPL_VIEW_CREATOR(unsigned long);
CPPYY_IMPL_VIEW_CREATOR(long long);
CPPYY_IMPL_VIEW_CREATOR(unsigned long long);
CPPYY_IMPL_VIEW_CREATOR(float);
CPPYY_IMPL_VIEW_CREATOR(double);
CPPYY_IMPL_VIEW_CREATOR(long double);
CPPYY_IMPL_VIEW_CREATOR(std::complex<float>);
CPPYY_IMPL_VIEW_CREATOR(std::complex<double>);
CPPYY_IMPL_VIEW_CREATOR(std::complex<int>);
CPPYY_IMPL_VIEW_CREATOR(std::complex<long>);

PyObject* CPyCppyy::CreateLowLevelView(const char** address, Py_ssize_t* shape) {
    return CreateLowLevelViewT<const char*>(address, shape);
}
