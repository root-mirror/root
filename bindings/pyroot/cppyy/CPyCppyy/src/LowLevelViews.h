#ifndef CPYCPPYY_LOWLEVELVIEWS_H
#define CPYCPPYY_LOWLEVELVIEWS_H

#include <complex>
#include <stddef.h>

#if __cplusplus > 201402L
#include <cstddef>
#endif


namespace CPyCppyy {

class Converter;

class LowLevelView {
public:
    PyObject_HEAD
    Py_buffer   fBufInfo;
    void**      fBuf;
    Converter*  fConverter;

public:
    void* get_buf() { return fBuf ? *fBuf : fBufInfo.buf; }
    void  set_buf(void** buf) { fBuf = buf; fBufInfo.buf = get_buf(); }
};

#define CPPYY_DECL_VIEW_CREATOR(type)                                        \
    PyObject* CreateLowLevelView(type*,  Py_ssize_t* shape = nullptr);       \
    PyObject* CreateLowLevelView(type**, Py_ssize_t* shape = nullptr)

CPPYY_DECL_VIEW_CREATOR(bool);
CPPYY_DECL_VIEW_CREATOR(signed char);
CPPYY_DECL_VIEW_CREATOR(unsigned char);
#if __cplusplus > 201402L
CPPYY_DECL_VIEW_CREATOR(std::byte);
#endif
CPPYY_DECL_VIEW_CREATOR(short);
CPPYY_DECL_VIEW_CREATOR(unsigned short);
CPPYY_DECL_VIEW_CREATOR(int);
CPPYY_DECL_VIEW_CREATOR(unsigned int);
CPPYY_DECL_VIEW_CREATOR(long);
CPPYY_DECL_VIEW_CREATOR(unsigned long);
CPPYY_DECL_VIEW_CREATOR(long long);
CPPYY_DECL_VIEW_CREATOR(unsigned long long);
CPPYY_DECL_VIEW_CREATOR(float);
CPPYY_DECL_VIEW_CREATOR(double);
CPPYY_DECL_VIEW_CREATOR(long double);
CPPYY_DECL_VIEW_CREATOR(std::complex<float>);
CPPYY_DECL_VIEW_CREATOR(std::complex<double>);
CPPYY_DECL_VIEW_CREATOR(std::complex<int>);
CPPYY_DECL_VIEW_CREATOR(std::complex<long>);

PyObject* CreateLowLevelView(const char**, Py_ssize_t* shape = nullptr);

inline PyObject* CreatePointerView(void* ptr, size_t size = (size_t)-1) {
    Py_ssize_t shape[] = {1, (Py_ssize_t)size};
    return CreateLowLevelView((uintptr_t*)ptr, shape);
}

//- low level view type and type verification --------------------------------
extern PyTypeObject LowLevelView_Type;

template<typename T>
inline bool LowLevelView_Check(T* object)
{
    return object && PyObject_TypeCheck(object, &LowLevelView_Type);
}

template<typename T>
inline bool LowLevelView_CheckExact(T* object)
{
    return object && Py_TYPE(object) == &LowLevelView_Type;
}

} // namespace CPyCppyy

#endif // !CPYCPPYY_LOWLEVELVIEWS_H
