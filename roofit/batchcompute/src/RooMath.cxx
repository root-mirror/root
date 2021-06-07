/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

// -- CLASS DESCRIPTION [MISC] --
// RooMath is a singleton class implementing various mathematical
// functions not found in TMath, mostly involving complex algebra
//
//
#include "RooMath.h"
#include "rbc.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>

namespace faddeeva_impl {
__device__ __host__ static inline void cexp(double& re, double& im)
    {
	// with gcc on unix machines and on x86_64, we can gain by hand-coding
	// exp(z) for the x87 coprocessor; other platforms have the default
	// routines as fallback implementation, and compilers other than gcc on
	// x86_64 generate better code with the default routines; also avoid
	// the inline assembly code when the copiler is not optimising code, or
	// is optimising for code size
	// (we insist on __unix__ here, since the assemblers on other OSs
	// typically do not speak AT&T syntax as gas does...)
#if defined(__CUDACC__) || !(defined(__GNUC__) || defined(__clang__)) || \
	!defined(__unix__) || !defined(__x86_64__) || \
	!defined(__OPTIMIZE__) || defined(__OPTIMIZE_SIZE__) || \
	defined(__INTEL_COMPILER) || \
	defined(__OPEN64__) || defined(__PATHSCALE__)
	const double e = std::exp(re);
	re = e * std::cos(im);
	im = e * std::sin(im);
#else
	__asm__ (
		"fxam\n\t"		 // examine st(0): NaN? Inf?
		"fstsw %%ax\n\t"
		"movb $0x45,%%dh\n\t"
		"andb %%ah,%%dh\n\t"
		"cmpb $0x05,%%dh\n\t"
		"jz 1f\n\t"		 // have NaN or infinity, handle below
		"fldl2e\n\t"		 // load log2(e)
		"fmulp\n\t"		 // re * log2(e)
		"fld %%st(0)\n\t"	 // duplicate re * log2(e)
		"frndint\n\t"		 // int(re * log2(e))
		"fsubr %%st,%%st(1)\n\t" // st(1) = x = frac(re * log2(e))
		"fxch\n\t"		 // swap st(0), st(1)
		"f2xm1\n\t"		 // 2^x - 1
		"fld1\n\t"		 // st(0) = 1
		"faddp\n\t"		 // st(0) = 2^x
		"fscale\n\t"		 // 2 ^ (int(re * log2(e)) + x)
		"fstp %%st(1)\n\t"	 // pop st(1)
		"jmp 2f\n\t"
		"1:\n\t"		 // handle NaN, Inf...
		"testl $0x200, %%eax\n\t"// -infinity?
		"jz 2f\n\t"
		"fstp %%st\n\t"		 // -Inf, so pop st(0)
		"fldz\n\t"		 // st(0) = 0
		"2:\n\t"		 // here. we have st(0) == exp(re)
		"fxch\n\t"		 // st(0) = im, st(1) = exp(re)
		"fsincos\n\t"		 // st(0) = cos(im), st(1) = sin(im)
		"fnstsw %%ax\n\t"
		"testl $0x400,%%eax\n\t"
		"jz 4f\n\t"		 // |im| too large for fsincos?
		"fldpi\n\t"		 // st(0) = pi
		"fadd %%st(0)\n\t"	 // st(0) *= 2;
		"fxch %%st(1)\n\t"	 // st(0) = im, st(1) = 2 * pi
		"3:\n\t"
		"fprem1\n\t"		 // st(0) = fmod(im, 2 * pi)
		"fnstsw %%ax\n\t"
		"testl $0x400,%%eax\n\t"
		"jnz 3b\n\t"		 // fmod done?
		"fstp %%st(1)\n\t"	 // yes, pop st(1) == 2 * pi
		"fsincos\n\t"		 // st(0) = cos(im), st(1) = sin(im)
		"4:\n\t"		 // all fine, fsincos succeeded
		"fmul %%st(2)\n\t"	 // st(0) *= st(2)
		"fxch %%st(2)\n\t"	 // st(2)=exp(re)*cos(im),st(0)=exp(im)
		"fmulp %%st(1)\n\t"	 // st(1)=exp(re)*sin(im), pop st(0)
		: "=t" (im), "=u" (re): "0" (re), "1" (im) :
		    "eax", "dh", "cc"
#ifndef __clang__
		    // normal compilers (like gcc) want to be told that we
		    // clobber x87 registers, even if we pop them afterwards
		    // (so they can make sure they don't save anything there)
		    , "st(5)", "st(6)", "st(7)"
#else // __clang__
		    // clang produces an error message with the clobber list
		    // from above - not sure why; it seems harmless to leave
		    // the popped x87 registers out of the clobber list for
		    // clang, and that is in fact what seems to be recommended
		    // here:
		    // http://lists.cs.uiuc.edu/pipermail/cfe-dev/2012-May/021715.html
#endif // __clang__
		    );
#endif
    }

    template <class T, unsigned N, unsigned NTAYLOR, unsigned NCF>
__device__ __host__ static inline std::complex<T> faddeeva_smabmq_impl(
	    T zre, T zim, const T tm,
	    const T (&a)[N], const T (&npi)[N],
	    const T (&taylorarr)[N * NTAYLOR * 2])
    {
	// catch singularities in the Fourier representation At
	// z = n pi / tm, and provide a Taylor series expansion in those
	// points, and only use it when we're close enough to the real axis
	// that there is a chance we need it
	const T zim2 = zim * zim;
	const T maxnorm = T(9) / T(1000000);
	if (zim2 < maxnorm) {
	    // we're close enough to the real axis that we need to worry about
	    // singularities
	    const T dnsing = tm * zre / npi[1];
	    const T dnsingmax2 = (T(N) - T(1) / T(2)) * (T(N) - T(1) / T(2));
	    if (dnsing * dnsing < dnsingmax2) {
		// we're in the interesting range of the real axis as well...
		// deal with Re(z) < 0 so we only need N different Taylor
		// expansions; use w(-x+iy) = conj(w(x+iy))
		const bool negrez = zre < T(0);
		// figure out closest singularity
		const int nsing = int(std::abs(dnsing) + T(1) / T(2));
		// and calculate just how far we are from it
		const T zmnpire = std::abs(zre) - npi[nsing];
		const T zmnpinorm = zmnpire * zmnpire + zim2;
		// close enough to one of the singularities?
		if (zmnpinorm < maxnorm) {
		    const T* coeffs = &taylorarr[nsing * NTAYLOR * 2];
		    // calculate value of taylor expansion...
		    // (note: there's no chance to vectorize this one, since
		    // the value of the next iteration depend on the ones from
		    // the previous iteration)
		    T sumre = coeffs[0], sumim = coeffs[1];
		    for (unsigned i = 1; i < NTAYLOR; ++i) {
			const T re = sumre * zmnpire - sumim * zim;
			const T im = sumim * zmnpire + sumre * zim;
			sumre = re + coeffs[2 * i + 0];
			sumim = im + coeffs[2 * i + 1];
		    }
		    // undo the flip in real part of z if needed
		    if (negrez) return std::complex<T>(sumre, -sumim);
		    else return std::complex<T>(sumre, sumim);
		}
	    }
	}
	// negative Im(z) is treated by calculating for -z, and using the
	// symmetry properties of erfc(z)
	const bool negimz = zim < T(0);
	if (negimz) {
	    zre = -zre;
	    zim = -zim;
	}
        const T znorm = zre * zre + zim2;
	if (znorm > tm * tm) {
	    // use continued fraction approximation for |z| large
	    const T isqrtpi = 5.64189583547756287e-01;
	    const T z2re = (zre + zim) * (zre - zim);
	    const T z2im = T(2) * zre * zim;
	    T cfre = T(1), cfim = T(0), cfnorm = T(1);
	    for (unsigned k = NCF; k; --k) {
		cfre = +(T(k) / T(2)) * cfre / cfnorm;
		cfim = -(T(k) / T(2)) * cfim / cfnorm;
		if (k & 1) cfre -= z2re, cfim -= z2im;
		else cfre += T(1);
		cfnorm = cfre * cfre + cfim * cfim;
	    }
	    T sumre =  (zim * cfre - zre * cfim) * isqrtpi / cfnorm;
	    T sumim = -(zre * cfre + zim * cfim) * isqrtpi / cfnorm;
	    if (negimz) {
		// use erfc(-z) = 2 - erfc(z) to get good accuracy for
		// Im(z) < 0: 2 / exp(z^2) - w(z)
		T ez2re = -z2re, ez2im = -z2im;
		faddeeva_impl::cexp(ez2re, ez2im);
		return std::complex<T>(T(2) * ez2re - sumre,
			T(2) * ez2im - sumim);
	    } else {
		return std::complex<T>(sumre, sumim);
	    }
	}
	const T twosqrtpi = 3.54490770181103205e+00;
	const T tmzre = tm * zre, tmzim = tm * zim;
	// calculate exp(i tm z)
	T eitmzre = -tmzim, eitmzim = tmzre;
	faddeeva_impl::cexp(eitmzre, eitmzim);
	// form 1 +/- exp (i tm z)
	const T numerarr[4] = {
	    T(1) - eitmzre, -eitmzim, T(1) + eitmzre, +eitmzim
	};
	// form tm z * (1 +/- exp(i tm z))
	const T numertmz[4] = {
	    tmzre * numerarr[0] - tmzim * numerarr[1],
	    tmzre * numerarr[1] + tmzim * numerarr[0],
	    tmzre * numerarr[2] - tmzim * numerarr[3],
	    tmzre * numerarr[3] + tmzim * numerarr[2]
	};
	// common subexpressions for use inside the loop
	const T reimtmzm2 = T(-2) * tmzre * tmzim;
	const T imtmz2 = tmzim * tmzim;
	const T reimtmzm22 = reimtmzm2 * reimtmzm2;
	// on non-x86_64 architectures, when the compiler is producing
	// unoptimised code and when optimising for code size, we use the
	// straightforward implementation, but for x86_64, we use the
	// brainf*cked code below that the gcc vectorizer likes to gain a few
	// clock cycles; non-gcc compilers also get the normal code, since they
	// usually do a better job with the default code (and yes, it's a pain
	// that they're all pretending to be gcc)
#if (!defined(__x86_64__)) || !defined(__OPTIMIZE__) || \
	defined(__OPTIMIZE_SIZE__) || defined(__INTEL_COMPILER) || \
	defined(__clang__) || defined(__OPEN64__) || \
	defined(__PATHSCALE__) || !defined(__GNUC__)
        T sumre = (-a[0] / znorm) * (numerarr[0] * zre + numerarr[1] * zim);
        T sumim = (-a[0] / znorm) * (numerarr[1] * zre - numerarr[0] * zim);
        for (unsigned i = 0; i < N; ++i) {
            const unsigned j = (i << 1) & 2;
            // denominator
            const T wk = imtmz2 + (npi[i] + tmzre) * (npi[i] - tmzre);
            // norm of denominator
            const T norm = wk * wk + reimtmzm22;
            const T f = T(2) * tm * a[i] / norm;
            // sum += a[i] * numer / wk
            sumre -= f * (numertmz[j] * wk + numertmz[j + 1] * reimtmzm2);
            sumim -= f * (numertmz[j + 1] * wk - numertmz[j] * reimtmzm2);
        }
#else
	// BEGIN fully vectorisable code - enjoy reading... ;)
	T tmp[2 * N];
	for (unsigned i = 0; i < N; ++i) {
	    const T wk = imtmz2 + (npi[i] + tmzre) * (npi[i] - tmzre);
	    tmp[2 * i + 0] = wk;
	    tmp[2 * i + 1] = T(2) * tm * a[i] / (wk * wk + reimtmzm22);
	}
	for (unsigned i = 0; i < N / 2; ++i) {
	    T wk = tmp[4 * i + 0], f = tmp[4 * i + 1];
	    tmp[4 * i + 0] = -f * (numertmz[0] * wk + numertmz[1] * reimtmzm2);
	    tmp[4 * i + 1] = -f * (numertmz[1] * wk - numertmz[0] * reimtmzm2);
	    wk = tmp[4 * i + 2], f = tmp[4 * i + 3];
	    tmp[4 * i + 2] = -f * (numertmz[2] * wk + numertmz[3] * reimtmzm2);
	    tmp[4 * i + 3] = -f * (numertmz[3] * wk - numertmz[2] * reimtmzm2);
	}
	if (N & 1) {
	    // we may have missed one element in the last loop; if so, process
	    // it now...
	    const T wk = tmp[2 * N - 2], f = tmp[2 * N - 1];
	    tmp[2 * (N - 1) + 0] = -f * (numertmz[0] * wk + numertmz[1] * reimtmzm2);
	    tmp[2 * (N - 1) + 1] = -f * (numertmz[1] * wk - numertmz[0] * reimtmzm2);
	}
	T sumre = (-a[0] / znorm) * (numerarr[0] * zre + numerarr[1] * zim);
	T sumim = (-a[0] / znorm) * (numerarr[1] * zre - numerarr[0] * zim);
	for (unsigned i = 0; i < N; ++i) {
	    sumre += tmp[2 * i + 0];
	    sumim += tmp[2 * i + 1];
	}
	// END fully vectorisable code
#endif
	// prepare the result
	if (negimz) {
	    // use erfc(-z) = 2 - erfc(z) to get good accuracy for
	    // Im(z) < 0: 2 / exp(z^2) - w(z)
	    const T z2im = -T(2) * zre * zim;
	    const T z2re = -(zre + zim) * (zre - zim);
	    T ez2re = z2re, ez2im = z2im;
	    faddeeva_impl::cexp(ez2re, ez2im);
	    return std::complex<T>(T(2) * ez2re + sumim / twosqrtpi,
		    T(2) * ez2im - sumre / twosqrtpi);
	} else {
	    return std::complex<T>(-sumim / twosqrtpi, sumre / twosqrtpi);
	}
    }

__device__ static constexpr double npi24[24] = { // precomputed values n * pi
	0.00000000000000000e+00, 3.14159265358979324e+00, 6.28318530717958648e+00,
	9.42477796076937972e+00, 1.25663706143591730e+01, 1.57079632679489662e+01,
	1.88495559215387594e+01, 2.19911485751285527e+01, 2.51327412287183459e+01,
	2.82743338823081391e+01, 3.14159265358979324e+01, 3.45575191894877256e+01,
	3.76991118430775189e+01, 4.08407044966673121e+01, 4.39822971502571053e+01,
	4.71238898038468986e+01, 5.02654824574366918e+01, 5.34070751110264851e+01,
	5.65486677646162783e+01, 5.96902604182060715e+01, 6.28318530717958648e+01,
	6.59734457253856580e+01, 6.91150383789754512e+01, 7.22566310325652445e+01,
    };
__device__ static constexpr double a24[24] = { // precomputed Fourier coefficient prefactors
	2.95408975150919338e-01, 2.75840233292177084e-01, 2.24573955224615866e-01,
	1.59414938273911723e-01, 9.86657664154541891e-02, 5.32441407876394120e-02,
	2.50521500053936484e-02, 1.02774656705395362e-02, 3.67616433284484706e-03,
	1.14649364124223317e-03, 3.11757015046197600e-04, 7.39143342960301488e-05,
	1.52794934280083635e-05, 2.75395660822107093e-06, 4.32785878190124505e-07,
	5.93003040874588103e-08, 7.08449030774820423e-09, 7.37952063581678038e-10,
	6.70217160600200763e-11, 5.30726516347079017e-12, 3.66432411346763916e-13,
	2.20589494494103134e-14, 1.15782686262855879e-15, 5.29871142946730482e-17,
    };
__device__ static constexpr double taylorarr24[24 * 12] = {
	// real part imaginary part, low order coefficients last
	// nsing = 0
	 0.00000000000000000e-00,  3.00901111225470020e-01,
	 5.00000000000000000e-01,  0.00000000000000000e-00,
	 0.00000000000000000e-00, -7.52252778063675049e-01,
	-1.00000000000000000e-00,  0.00000000000000000e-00,
	 0.00000000000000000e-00,  1.12837916709551257e+00,
	 1.00000000000000000e-00,  0.00000000000000000e-00,
	// nsing = 1
	-2.22423508493755319e-01,  1.87966717746229718e-01,
	 3.41805419240637628e-01,  3.42752593807919263e-01,
	 4.66574321730757753e-01, -5.59649213591058097e-01,
	-8.05759710273191021e-01, -5.38989366115424093e-01,
	-4.88914083733395200e-01,  9.80580906465856792e-01,
	 9.33757118080975970e-01,  2.82273885115127769e-01,
	// nsing = 2
	-2.60522586513312894e-01, -4.26259455096092786e-02,
	 1.36549702008863349e-03,  4.39243227763478846e-01,
	 6.50591493715480700e-01, -1.23422352472779046e-01,
	-3.43379903564271318e-01, -8.13862662890748911e-01,
	-7.96093943501906645e-01,  6.11271022503935772e-01,
	 7.60213717643090957e-01,  4.93801903948967945e-01,
	// nsing = 3
	-1.18249853727020186e-01, -1.90471659765411376e-01,
	-2.59044664869706839e-01,  2.69333898502392004e-01,
	 4.99077838344125714e-01,  2.64644800189075006e-01,
	 1.26114512111568737e-01, -7.46519337025968199e-01,
	-8.47666863706379907e-01,  1.89347715957263646e-01,
	 5.39641485816297176e-01,  5.97805988669631615e-01,
	// nsing = 4
	 4.94825297066481491e-02, -1.71428212158876197e-01,
	-2.97766677111471585e-01,  1.60773286596649656e-02,
	 1.88114210832460682e-01,  4.11734391195006462e-01,
	 3.98540613293909842e-01, -4.63321903522162715e-01,
	-6.99522070542463639e-01, -1.32412024008354582e-01,
	 3.33997185986131785e-01,  6.01983450812696742e-01,
	// nsing = 5
	 1.18367078448232332e-01, -6.09533063579086850e-02,
	-1.74762998833038991e-01, -1.39098099222000187e-01,
	-6.71534655984154549e-02,  3.34462251996496680e-01,
	 4.37429678577360024e-01, -1.59613865629038012e-01,
	-4.71863911886034656e-01, -2.92759316465055762e-01,
	 1.80238737704018306e-01,  5.42834914744283253e-01,
	// nsing = 6
	 8.87698096005701290e-02,  2.84339354980994902e-02,
	-3.18943083830766399e-02, -1.53946887977045862e-01,
	-1.71825061547624858e-01,  1.70734367410600348e-01,
	 3.33690792296469441e-01,  3.97048587678703930e-02,
	-2.66422678503135697e-01, -3.18469797424381480e-01,
	 8.48049724711137773e-02,  4.60546329221462864e-01,
	// nsing = 7
	 2.99767046276705077e-02,  5.34659695701718247e-02,
	 4.53131030251822568e-02, -9.37915401977138648e-02,
	-1.57982359988083777e-01,  3.82170507060760740e-02,
	 1.98891589845251706e-01,  1.17546677047049354e-01,
	-1.27514335237079297e-01, -2.72741112680307074e-01,
	 3.47906344595283767e-02,  3.82277517244493224e-01,
	// nsing = 8
	-7.35922494437203395e-03,  3.72011290318534610e-02,
	 5.66783220847204687e-02, -3.21015398169199501e-02,
	-1.00308737825172555e-01, -2.57695148077963515e-02,
	 9.67294850588435368e-02,  1.18174625238337507e-01,
	-5.21266530264988508e-02, -2.08850084114630861e-01,
	 1.24443217440050976e-02,  3.19239968065752286e-01,
	// nsing = 9
	-1.66126772808035320e-02,  1.46180329587665321e-02,
	 3.85927576915247303e-02,  1.18910471133003227e-03,
	-4.94003498320899806e-02, -3.93468443660139110e-02,
	 3.92113167048952835e-02,  9.03306084789976219e-02,
	-1.82889636251263500e-02, -1.53816215444915245e-01,
	 3.88103861995563741e-03,  2.72090310854550347e-01,
	// nsing = 10
	-1.21245068916826880e-02,  1.59080224420074489e-03,
	 1.91116222508366035e-02,  1.05879549199053302e-02,
	-1.97228428219695318e-02, -3.16962067712639397e-02,
	 1.34110372628315158e-02,  6.18045654429108837e-02,
	-5.52574921865441838e-03, -1.14259663804569455e-01,
	 1.05534036292203489e-03,  2.37326534898818288e-01,
	// nsing = 11
	-5.96835002183177493e-03, -2.42594931567031205e-03,
	 7.44753817476594184e-03,  9.33450807578394386e-03,
	-6.52649522783026481e-03, -2.08165802069352019e-02,
	 3.89988065678848650e-03,  4.12784313451549132e-02,
	-1.44110721106127920e-03, -8.76484782997757425e-02,
	 2.50210184908121337e-04,  2.11131066219336647e-01,
	// nsing = 12
	-2.24505212235034193e-03, -2.38114524227619446e-03,
	 2.36375918970809340e-03,  5.97324040603806266e-03,
	-1.81333819936645381e-03, -1.28126250720444051e-02,
	 9.69251586187208358e-04,  2.83055679874589732e-02,
	-3.24986363596307374e-04, -6.97056268370209313e-02,
	 5.17231862038123061e-05,  1.90681117197597520e-01,
	// nsing = 13
	-6.76887607549779069e-04, -1.48589685249767064e-03,
	 6.22548369472046953e-04,  3.43871156746448680e-03,
	-4.26557147166379929e-04, -7.98854145009655400e-03,
	 2.06644460919535524e-04,  2.03107152586353217e-02,
	-6.34563929410856987e-05, -5.71425144910115832e-02,
	 9.32252179140502456e-06,  1.74167663785025829e-01,
	// nsing = 14
	-1.67596437777156162e-04, -8.05384193869903178e-04,
	 1.37627277777023791e-04,  1.97652692602724093e-03,
	-8.54392244879459717e-05, -5.23088906415977167e-03,
	 3.78965577556493513e-05,  1.52191559129376333e-02,
	-1.07393019498185646e-05, -4.79347862153366295e-02,
	 1.46503970628861795e-06,  1.60471011683477685e-01,
	// nsing = 15
	-3.45715760630978778e-05, -4.31089554210205493e-04,
	 2.57350138106549737e-05,  1.19449262097417514e-03,
	-1.46322227517372253e-05, -3.61303766799909378e-03,
	 5.99057675687392260e-06,  1.17993805017130890e-02,
	-1.57660578509526722e-06, -4.09165023743669707e-02,
	 2.00739683204152177e-07,  1.48879348585662670e-01,
	// nsing = 16
	-5.99735188857573424e-06, -2.42949218855805052e-04,
	 4.09249090936269722e-06,  7.67400152727128171e-04,
	-2.14920611287648034e-06, -2.60710519575546230e-03,
	 8.17591694958640978e-07,  9.38581640137393053e-03,
	-2.00910914042737743e-07, -3.54045580123653803e-02,
	 2.39819738182594508e-08,  1.38916449405613711e-01,
	// nsing = 17
	-8.80708505155966658e-07, -1.46479474515521504e-04,
	 5.55693207391871904e-07,  5.19165587844615415e-04,
	-2.71391142598826750e-07, -1.94439427580099576e-03,
	 9.64641799864928425e-08,  7.61536975207357980e-03,
	-2.22357616069432967e-08, -3.09762939485679078e-02,
	 2.49806920458212581e-09,  1.30247401712293206e-01,
	// nsing = 18
	-1.10007111030476390e-07, -9.35886150886691786e-05,
	 6.46244096997824390e-08,  3.65267193418479043e-04,
	-2.95175785569292542e-08, -1.48730955943961081e-03,
	 9.84949251974795537e-09,  6.27824679148707177e-03,
	-2.13827217704781576e-09, -2.73545766571797965e-02,
	 2.26877724435352177e-10,  1.22627158810895267e-01,
	// nsing = 19
	-1.17302439957657553e-08, -6.24890956722053332e-05,
	 6.45231881609786173e-09,  2.64799907072561543e-04,
	-2.76943921343331654e-09, -1.16094187847598385e-03,
	 8.71074689656480749e-10,  5.24514377390761210e-03,
	-1.78730768958639407e-10, -2.43489203319091538e-02,
	 1.79658223341365988e-11,  1.15870972518909888e-01,
	// nsing = 20
	-1.07084502471985403e-09, -4.31515421260633319e-05,
	 5.54152563270547927e-10,  1.96606443937168357e-04,
	-2.24423474431542338e-10, -9.21550077887211094e-04,
	 6.67734377376211580e-11,  4.43201203646827019e-03,
	-1.29896907717633162e-11, -2.18236356404862774e-02,
	 1.24042409733678516e-12,  1.09836276968151848e-01,
	// nsing = 21
	-8.38816525569060600e-11, -3.06091807093959821e-05,
	 4.10033961556230842e-11,  1.48895624771753491e-04,
	-1.57238128435253905e-11, -7.42073499862065649e-04,
	 4.43938379112418832e-12,  3.78197089773957382e-03,
	-8.21067867869285873e-13, -1.96793607299577220e-02,
	 7.46725770201828754e-14,  1.04410965521273064e-01,
	// nsing = 22
	-5.64848922712870507e-12, -2.22021942382507691e-05,
	 2.61729537775838587e-12,  1.14683068921649992e-04,
	-9.53316139085394895e-13, -6.05021573565916914e-04,
	 2.56116039498542220e-13,  3.25530796858307225e-03,
	-4.51482829896525004e-14, -1.78416955716514289e-02,
	 3.91940313268087086e-15,  9.95054815464739996e-02,
	// nsing = 23
	-3.27482357793897640e-13, -1.64138890390689871e-05,
	 1.44278798346454523e-13,  8.96362542918265398e-05,
	-5.00524303437266481e-14, -4.98699756861136127e-04,
	 1.28274026095767213e-14,  2.82359118537843949e-03,
	-2.16009593993917109e-15, -1.62538825704327487e-02,
	 1.79368667683853708e-16,  9.50473084594884184e-02
    };

    const double npi11[11] = { // precomputed values n * pi
	0.00000000000000000e+00, 3.14159265358979324e+00, 6.28318530717958648e+00,
	9.42477796076937972e+00, 1.25663706143591730e+01, 1.57079632679489662e+01,
	1.88495559215387594e+01, 2.19911485751285527e+01, 2.51327412287183459e+01,
	2.82743338823081391e+01, 3.14159265358979324e+01
    };
    const double a11[11] = { // precomputed Fourier coefficient prefactors
	4.43113462726379007e-01, 3.79788034073635143e-01, 2.39122407410867584e-01,
	1.10599187402169792e-01, 3.75782250080904725e-02, 9.37936104296856288e-03,
	1.71974046186334976e-03, 2.31635559000523461e-04, 2.29192401420125452e-05,
	1.66589592139340077e-06, 8.89504561311882155e-08
    };
    const double taylorarr11[11 * 6] = {
	// real part imaginary part, low order coefficients last
	// nsing = 0
	-1.00000000000000000e+00,  0.00000000000000000e+00,
	 0.00000000000000000e-01,  1.12837916709551257e+00,
	 1.00000000000000000e+00,  0.00000000000000000e+00,
	// nsing = 1
	-5.92741768247463996e-01, -7.19914991991294310e-01,
	-6.73156763521649944e-01,  8.14025039279059577e-01,
	 8.57089811121701143e-01,  4.00248106586639754e-01,
	// nsing = 2
	 1.26114512111568737e-01, -7.46519337025968199e-01,
	-8.47666863706379907e-01,  1.89347715957263646e-01,
	 5.39641485816297176e-01,  5.97805988669631615e-01,
	// nsing = 3
	 4.43238482668529408e-01, -3.03563167310638372e-01,
	-5.88095866853990048e-01, -2.32638360700858412e-01,
	 2.49595637924601714e-01,  5.77633779156009340e-01,
	// nsing = 4
	 3.33690792296469441e-01,  3.97048587678703930e-02,
	-2.66422678503135697e-01, -3.18469797424381480e-01,
	 8.48049724711137773e-02,  4.60546329221462864e-01,
	// nsing = 5
	 1.42043544696751869e-01,  1.24094227867032671e-01,
	-8.31224229982140323e-02, -2.40766729258442100e-01,
	 2.11669512031059302e-02,  3.48650139549945097e-01,
	// nsing = 6
	 3.92113167048952835e-02,  9.03306084789976219e-02,
	-1.82889636251263500e-02, -1.53816215444915245e-01,
	 3.88103861995563741e-03,  2.72090310854550347e-01,
	// nsing = 7
	 7.37741897722738503e-03,  5.04625223970221539e-02,
	-2.87394336989990770e-03, -9.96122819257496929e-02,
	 5.22745478269428248e-04,  2.23361039070072101e-01,
	// nsing = 8
	 9.69251586187208358e-04,  2.83055679874589732e-02,
	-3.24986363596307374e-04, -6.97056268370209313e-02,
	 5.17231862038123061e-05,  1.90681117197597520e-01,
	// nsing = 9
	 9.01625563468897100e-05,  1.74961124275657019e-02,
	-2.65745127697337342e-05, -5.22070356354932341e-02,
	 3.75952450449939411e-06,  1.67018782142871146e-01,
	// nsing = 10
	 5.99057675687392260e-06,  1.17993805017130890e-02,
	-1.57660578509526722e-06, -4.09165023743669707e-02,
	 2.00739683204152177e-07,  1.48879348585662670e-01
    };
}

__device__ __host__ std::complex<double> RooMath::faddeeva(std::complex<double> z)
{
    return faddeeva_impl::faddeeva_smabmq_impl<double, 24, 6, 9>(
	    z.real(), z.imag(), 12., faddeeva_impl::a24,
	    faddeeva_impl::npi24, faddeeva_impl::taylorarr24);
}

std::complex<double> RooMath::faddeeva_fast(std::complex<double> z)
{
    return faddeeva_impl::faddeeva_smabmq_impl<double, 11, 3, 3>(
	    z.real(), z.imag(), 8., faddeeva_impl::a11,
	    faddeeva_impl::npi11, faddeeva_impl::taylorarr11);
}

std::complex<double> RooMath::erfc(const std::complex<double> z)
{
    double re = -z.real() * z.real() + z.imag() * z.imag();
    double im = -2. * z.real() * z.imag();
    faddeeva_impl::cexp(re, im);
    return (z.real() >= 0.) ?
	(std::complex<double>(re, im) *
	 faddeeva(std::complex<double>(-z.imag(), z.real()))) :
	(2. - std::complex<double>(re, im) *
	 faddeeva(std::complex<double>(z.imag(), -z.real())));
}

std::complex<double> RooMath::erfc_fast(const std::complex<double> z)
{
    double re = -z.real() * z.real() + z.imag() * z.imag();
    double im = -2. * z.real() * z.imag();
    faddeeva_impl::cexp(re, im);
    return (z.real() >= 0.) ?
	(std::complex<double>(re, im) *
	 faddeeva_fast(std::complex<double>(-z.imag(), z.real()))) :
	(2. - std::complex<double>(re, im) *
	 faddeeva_fast(std::complex<double>(z.imag(), -z.real())));
}

std::complex<double> RooMath::erf(const std::complex<double> z)
{
    double re = -z.real() * z.real() + z.imag() * z.imag();
    double im = -2. * z.real() * z.imag();
    faddeeva_impl::cexp(re, im);
    return (z.real() >= 0.) ?
	(1. - std::complex<double>(re, im) *
	 faddeeva(std::complex<double>(-z.imag(), z.real()))) :
	(std::complex<double>(re, im) *
	 faddeeva(std::complex<double>(z.imag(), -z.real())) - 1.);
}

std::complex<double> RooMath::erf_fast(const std::complex<double> z)
{
    double re = -z.real() * z.real() + z.imag() * z.imag();
    double im = -2. * z.real() * z.imag();
    faddeeva_impl::cexp(re, im);
    return (z.real() >= 0.) ?
	(1. - std::complex<double>(re, im) *
	 faddeeva_fast(std::complex<double>(-z.imag(), z.real()))) :
	(std::complex<double>(re, im) *
	 faddeeva_fast(std::complex<double>(z.imag(), -z.real())) - 1.);
}


Double_t RooMath::interpolate(Double_t ya[], Int_t n, Double_t x) 
{
  // Interpolate array 'ya' with 'n' elements for 'x' (between 0 and 'n'-1)

  // Int to Double conversion is faster via array lookup than type conversion!
  static Double_t itod[20] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
			      10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0} ;
  int i,m,ns=1 ;
  Double_t den,dif,dift/*,ho,hp,w*/,y,dy ;
  Double_t c[20], d[20] ;

  dif = fabs(x) ;
  for(i=1 ; i<=n ; i++) {
    if ((dift=fabs(x-itod[i-1]))<dif) {
      ns=i ;
      dif=dift ;
    }
    c[i] = ya[i-1] ;
    d[i] = ya[i-1] ;
  }
  
  y=ya[--ns] ;
  for(m=1 ; m<n; m++) {       
    for(i=1 ; i<=n-m ; i++) { 
      den=(c[i+1] - d[i])/itod[m] ;
      d[i]=(x-itod[i+m-1])*den ;
      c[i]=(x-itod[i-1])*den ;
    }
    dy = (2*ns)<(n-m) ? c[ns+1] : d[ns--] ;
    y += dy ;
  }
  return y ;
}



Double_t RooMath::interpolate(Double_t xa[], Double_t ya[], Int_t n, Double_t x) 
{
  // Interpolate array 'ya' with 'n' elements for 'xa' 

  // Int to Double conversion is faster via array lookup than type conversion!
//   static Double_t itod[20] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
// 			      10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0} ;
  int i,m,ns=1 ;
  Double_t den,dif,dift,ho,hp,w,y,dy ;
  Double_t c[20], d[20] ;

  dif = fabs(x-xa[0]) ;
  for(i=1 ; i<=n ; i++) {
    if ((dift=fabs(x-xa[i-1]))<dif) {
      ns=i ;
      dif=dift ;
    }
    c[i] = ya[i-1] ;
    d[i] = ya[i-1] ;
  }
  
  y=ya[--ns] ;
  for(m=1 ; m<n; m++) {       
    for(i=1 ; i<=n-m ; i++) { 
      ho=xa[i-1]-x ;
      hp=xa[i-1+m]-x ;
      w=c[i+1]-d[i] ;
      den=ho-hp ;
      if (den==0.) {
        std::cerr << "In " << __func__ << "(), " << __FILE__ << ":" << __LINE__ << 
          ": Zero distance between points not allowed." << std::endl;
        return 0 ;
      }
      den = w/den ;
      d[i]=hp*den ;
      c[i]=ho*den;
    }
    dy = (2*ns)<(n-m) ? c[ns+1] : d[ns--] ;
    y += dy ;
  }
  return y ;
}
