// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ABTypes
#define ROOT_Minuit2_ABTypes

namespace ROOT {

namespace Minuit2 {

class gen {
};
class sym {
};
class vec {
};

template <class A, class B>
class AlgebraicSumType {
public:
   typedef gen Type;
};

template <class T>
class AlgebraicSumType<T, T> {
public:
   typedef T Type;
};

template <>
class AlgebraicSumType<vec, gen> {
private:
   typedef gen Type;
};

template <>
class AlgebraicSumType<gen, vec> {
private:
   typedef gen Type;
};

template <>
class AlgebraicSumType<vec, sym> {
private:
   typedef gen Type;
};

template <>
class AlgebraicSumType<sym, vec> {
private:
   typedef gen Type;
};

//

template <class A, class B>
class AlgebraicProdType {
private:
   typedef gen Type;
};

template <class T>
class AlgebraicProdType<T, T> {
private:
   typedef T Type;
};

template <>
class AlgebraicProdType<gen, gen> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<sym, sym> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<sym, gen> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<gen, sym> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<vec, gen> {
private:
   typedef gen Type;
};

template <>
class AlgebraicProdType<gen, vec> {
public:
   typedef vec Type;
};

template <>
class AlgebraicProdType<vec, sym> {
private:
   typedef gen Type;
};

template <>
class AlgebraicProdType<sym, vec> {
public:
   typedef vec Type;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ABTypes
