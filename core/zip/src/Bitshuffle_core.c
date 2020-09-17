/*
 * Bitshuffle - Filter for improving compression of typed binary data.
 *
 * Author: Kiyoshi Masui <kiyo@physics.ubc.ca>
 * Website: http://www.github.com/kiyo-masui/bitshuffle
 * Created: 2014
 *
 * See LICENSE file for details about copyright and rights to use.
 *
 */

#include "Bitshuffle_core.h"
#include "Bitshuffle_internals.h"

#include <stdio.h>
#include <string.h>


// Conditional includes for SSE2 and AVX2.
#if defined(_MSC_VER)
     /* Microsoft C/C++-compatible compiler */
     #include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
     /* GCC-compatible compiler, targeting x86/x86-64 */
     #include <x86intrin.h>
#elif defined(__GNUC__) && defined(__ARM_NEON__)
     /* GCC-compatible compiler, targeting ARM with NEON */
     #include <arm_neon.h>
#elif defined(__GNUC__) && defined(__IWMMXT__)
     /* GCC-compatible compiler, targeting ARM with WMMX */
     #include <mmintrin.h>
#elif (defined(__GNUC__) || defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
     /* XLC or GCC-compatible compiler, targeting PowerPC with VMX/VSX */
     #include <altivec.h>
#elif defined(__GNUC__) && defined(__SPE__)
     /* GCC-compatible compiler, targeting PowerPC with SPE */
     #include <spe.h>
#elif defined(__clang__) && defined(__linux__)
    #include <x86intrin.h>
#endif

#if defined (__x86_64__) && defined (__linux__)
#include "cpuid.h"
#endif

// from cpuid.h
#ifndef bit_AVX2
# define bit_AVX2 0x00000020
#endif

#ifndef bit_SSE4_2
# define bit_SSE4_2 0x100000
#endif

#ifndef bit_SSE2
# define bit_SSE2 0x4000000
#endif

#if defined(_MSC_VER)
typedef int64_t omp_size_t;
#else
typedef size_t omp_size_t;
#endif

// Macros.
#define CHECK_MULT_EIGHT(n) if ((n) % 8) return -80;
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


/* ---- Worker code not requiring special instruction sets. ----
 *
 * The following code does not use any x86 specific vectorized instructions
 * and should compile on any machine
 *
 */

/* Transpose 8x8 bit array packed into a single quadword *x*.
 * *t* is workspace. */
#define TRANS_BIT_8X8(x, t) {                                               \
        (t) = ((x) ^ ((x) >> 7)) & 0x00AA00AA00AA00AALL;                          \
        (x) = (x) ^ (t) ^ ((t) << 7);                                               \
        (t) = ((x) ^ ((x) >> 14)) & 0x0000CCCC0000CCCCLL;                         \
        (x) = (x) ^ (t) ^ ((t) << 14);                                              \
        (t) = ((x) ^ ((x) >> 28)) & 0x00000000F0F0F0F0LL;                         \
        (x) = (x) ^ (t) ^ ((t) << 28);                                              \
    }

/* Transpose 8x8 bit array along the diagonal from upper right
   to lower left */
#define TRANS_BIT_8X8_BE(x, t) {                                            \
        (t) = ((x) ^ ((x) >> 9)) & 0x0055005500550055LL;                          \
        (x) = (x) ^ (t) ^ ((t) << 9);                                               \
        (t) = ((x) ^ ((x) >> 18)) & 0x0000333300003333LL;                         \
        (x) = (x) ^ (t) ^ ((t) << 18);                                              \
        (t) = ((x) ^ ((x) >> 36)) & 0x000000000F0F0F0FLL;                         \
        (x) = (x) ^ (t) ^ ((t) << 36);                                              \
    }

/* Transpose of an array of arbitrarily typed elements. */
#define TRANS_ELEM_TYPE(in, out, lda, ldb, type_t) {                        \
        size_t ii, jj, kk;                                                  \
        const type_t* in_type = (const type_t*) (in);                                 \
        type_t* out_type = (type_t*) (out);                                   \
        for(ii = 0; ii + 7 < (lda); ii += 8) {                                \
            for(jj = 0; jj < (ldb); jj++) {                                   \
                for(kk = 0; kk < 8; kk++) {                                 \
                    out_type[jj*(lda) + ii + kk] =                            \
                        in_type[ii*(ldb) + kk * (ldb) + jj];                    \
                }                                                           \
            }                                                               \
        }                                                                   \
        for(ii = (lda) - (lda) % 8; ii < (lda); ii ++) {                          \
            for(jj = 0; jj < (ldb); jj++) {                                   \
                out_type[jj*(lda) + ii] = in_type[ii*(ldb) + jj];                            \
            }                                                               \
        }                                                                   \
    }


/* Memory copy with bshuf call signature. For testing and profiling. */
int64_t bshuf_copy(const void* in, void* out, const size_t size,
         const size_t elem_size) {

    const char* in_b = (const char*) in;
    char* out_b = (char*) out;

    memcpy(out_b, in_b, size * elem_size);
    return size * elem_size;
}


/* Transpose bytes within elements, starting partway through input. */
int64_t bshuf_trans_byte_elem_remainder(const void* in, void* out, const size_t size,
         const size_t elem_size, const size_t start) {

    size_t ii, jj, kk;
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;

    CHECK_MULT_EIGHT(start);

    if (size > start) {
        // ii loop separated into 2 loops so the compiler can unroll
        // the inner one.
        for (ii = start; ii + 7 < size; ii += 8) {
            for (jj = 0; jj < elem_size; jj++) {
                for (kk = 0; kk < 8; kk++) {
                    out_b[jj * size + ii + kk]
                        = in_b[ii * elem_size + kk * elem_size + jj];
                }
            }
        }
        for (ii = size - size % 8; ii < size; ii ++) {
            for (jj = 0; jj < elem_size; jj++) {
                out_b[jj * size + ii] = in_b[ii * elem_size + jj];
            }
        }
    }
    return size * elem_size;
}


/* Transpose bytes within elements. */
int64_t bshuf_trans_byte_elem_scal(const void* in, void* out, const size_t size,
         const size_t elem_size) {

    return bshuf_trans_byte_elem_remainder(in, out, size, elem_size, 0);
}


/* Transpose bits within bytes. */
int64_t bshuf_trans_bit_byte_remainder(const void* in, void* out, const size_t size,
         const size_t elem_size, const size_t start_byte) {

    const uint64_t* in_b = (const uint64_t*) in;
    uint8_t* out_b = (uint8_t*) out;

    uint64_t x, t;

    size_t ii, kk;
    size_t nbyte = elem_size * size;
    size_t nbyte_bitrow = nbyte / 8;

    uint64_t e=1;
    const int little_endian = *(uint8_t *) &e == 1;
    const size_t bit_row_skip = little_endian ? nbyte_bitrow : -nbyte_bitrow;
    const int64_t bit_row_offset = little_endian ? 0 : 7 * nbyte_bitrow;

    CHECK_MULT_EIGHT(nbyte);
    CHECK_MULT_EIGHT(start_byte);

    for (ii = start_byte / 8; ii < nbyte_bitrow; ii ++) {
        x = in_b[ii];
        if (little_endian) {
            TRANS_BIT_8X8(x, t);
        } else {
            TRANS_BIT_8X8_BE(x, t);
        }
        for (kk = 0; kk < 8; kk ++) {
            out_b[bit_row_offset + kk * bit_row_skip + ii] = x;
            x = x >> 8;
        }
    }
    return size * elem_size;
}


/* Transpose bits within bytes. */
int64_t bshuf_trans_bit_byte_scal(const void* in, void* out, const size_t size,
         const size_t elem_size) {

    return bshuf_trans_bit_byte_remainder(in, out, size, elem_size, 0);
}


/* General transpose of an array, optimized for large element sizes. */
int64_t bshuf_trans_elem(const void* in, void* out, const size_t lda,
        const size_t ldb, const size_t elem_size) {

    size_t ii, jj;
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;
    for(ii = 0; ii < lda; ii++) {
        for(jj = 0; jj < ldb; jj++) {
            memcpy(&out_b[(jj*lda + ii) * elem_size],
                   &in_b[(ii*ldb + jj) * elem_size], elem_size);
        }
    }
    return lda * ldb * elem_size;
}


/* Transpose rows of shuffled bits (size / 8 bytes) within groups of 8. */
int64_t bshuf_trans_bitrow_eight(const void* in, void* out, const size_t size,
         const size_t elem_size) {

    size_t nbyte_bitrow = size / 8;

    CHECK_MULT_EIGHT(size);

    return bshuf_trans_elem(in, out, 8, elem_size, nbyte_bitrow);
}


/* Transpose bits within elements. */
int64_t bshuf_trans_bit_elem_scal(const void* in, void* out, const size_t size,
         const size_t elem_size) {

    int64_t count;
    void *tmp_buf;

    CHECK_MULT_EIGHT(size);

    tmp_buf = malloc(size * elem_size);
    if (tmp_buf == NULL) return -1;

    count = bshuf_trans_byte_elem_scal(in, out, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count = bshuf_trans_bit_byte_scal(out, tmp_buf, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count = bshuf_trans_bitrow_eight(tmp_buf, out, size, elem_size);

    free(tmp_buf);

    return count;
}


/* For data organized into a row for each bit (8 * elem_size rows), transpose
 * the bytes. */
int64_t bshuf_trans_byte_bitrow_scal(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    size_t ii, jj, kk, nbyte_row;
    const char *in_b;
    char *out_b;


    in_b = (const char*) in;
    out_b = (char*) out;

    nbyte_row = size / 8;

    CHECK_MULT_EIGHT(size);

    for (jj = 0; jj < elem_size; jj++) {
        for (ii = 0; ii < nbyte_row; ii++) {
            for (kk = 0; kk < 8; kk++) {
                out_b[ii * 8 * elem_size + jj * 8 + kk] = \
                        in_b[(jj * 8 + kk) * nbyte_row + ii];
            }
        }
    }
    return size * elem_size;
}


/* Shuffle bits within the bytes of eight element blocks. */
int64_t bshuf_shuffle_bit_eightelem_scal(const void* in, void* out, \
        const size_t size, const size_t elem_size) {

    const char *in_b;
    char *out_b;
    uint64_t x, t;
    size_t ii, jj, kk;
    size_t nbyte, out_index;

    uint64_t e=1;
    const int little_endian = *(uint8_t *) &e == 1;
    const size_t elem_skip = little_endian ? elem_size : -elem_size;
    const uint64_t elem_offset = little_endian ? 0 : 7 * elem_size;

    CHECK_MULT_EIGHT(size);

    in_b = (const char*) in;
    out_b = (char*) out;

    nbyte = elem_size * size;

    for (jj = 0; jj + 7 < 8 * elem_size; jj += 8) {
        for (ii = 0; ii + 8 * elem_size - 1 < nbyte; ii += 8 * elem_size) {
            x = *((uint64_t*) &in_b[ii + jj]);
            if (little_endian) {
                TRANS_BIT_8X8(x, t);
            } else {
                TRANS_BIT_8X8_BE(x, t);
            }
            for (kk = 0; kk < 8; kk++) {
                out_index = ii + jj / 8 + elem_offset + kk * elem_skip;
                *((uint8_t*) &out_b[out_index]) = x;
                x = x >> 8;
            }
        }
    }
    return size * elem_size;
}


/* Untranspose bits within elements. */
int64_t bshuf_untrans_bit_elem_scal(const void* in, void* out, const size_t size,
         const size_t elem_size) {

    int64_t count;
    void *tmp_buf;

    CHECK_MULT_EIGHT(size);

    tmp_buf = malloc(size * elem_size);
    if (tmp_buf == NULL) return -1;

    count = bshuf_trans_byte_bitrow_scal(in, tmp_buf, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count =  bshuf_shuffle_bit_eightelem_scal(tmp_buf, out, size, elem_size);

    free(tmp_buf);

    return count;
}

#if defined (__x86_64__) && defined (__linux__) && ((__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 9)) || (__clang__))

#ifdef _MSC_VER

/* MSC doesn't have __builtin_expect.  Just ignore likely/unlikely and
 * hope the compiler optimizes for the best. */

#define likely(x)       (x)
#define unlikely(x)     (x)

#else

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

#endif


/* ---- Worker code that uses SSE2 ----
 *
 * The following code makes use of the SSE2 instruction set and specialized
 * 16 byte registers. The SSE2 instructions are present on modern x86
 * processors. The first Intel processor microarchitecture supporting SSE2 was
 * Pentium 4 (2000).
 *
 */


/* Transpose bytes within elements for 16 bit elements. */
__attribute__ ((target ("sse2")))
int64_t bshuf_trans_byte_elem_SSE_16(const void* in, void* out, const size_t size) {
    
    size_t ii;
    const char *in_b = (const char*) in;
    char *out_b = (char*) out;
    __m128i a0, b0, a1, b1;

    for (ii=0; ii + 15 < size; ii += 16) {
        a0 = _mm_loadu_si128((__m128i *) &in_b[2*ii + 0*16]);
        b0 = _mm_loadu_si128((__m128i *) &in_b[2*ii + 1*16]);

        a1 = _mm_unpacklo_epi8(a0, b0);
        b1 = _mm_unpackhi_epi8(a0, b0);

        a0 = _mm_unpacklo_epi8(a1, b1);
        b0 = _mm_unpackhi_epi8(a1, b1);

        a1 = _mm_unpacklo_epi8(a0, b0);
        b1 = _mm_unpackhi_epi8(a0, b0);

        a0 = _mm_unpacklo_epi8(a1, b1);
        b0 = _mm_unpackhi_epi8(a1, b1);

        _mm_storeu_si128((__m128i *) &out_b[0*size + ii], a0);
        _mm_storeu_si128((__m128i *) &out_b[1*size + ii], b0);
    }
    return bshuf_trans_byte_elem_remainder(in, out, size, 2,
            size - size % 16);
}


/* Transpose bytes within elements for 32 bit elements. */
__attribute__ ((target ("sse2")))
int64_t bshuf_trans_byte_elem_SSE_32(const void* in, void* out, const size_t size) {
    
    size_t ii;
    const char *in_b;
    char *out_b;
    in_b = (const char*) in;
    out_b = (char*) out;
    __m128i a0, b0, c0, d0, a1, b1, c1, d1;

    for (ii=0; ii + 15 < size; ii += 16) {
        a0 = _mm_loadu_si128((__m128i *) &in_b[4*ii + 0*16]);
        b0 = _mm_loadu_si128((__m128i *) &in_b[4*ii + 1*16]);
        c0 = _mm_loadu_si128((__m128i *) &in_b[4*ii + 2*16]);
        d0 = _mm_loadu_si128((__m128i *) &in_b[4*ii + 3*16]);

        a1 = _mm_unpacklo_epi8(a0, b0);
        b1 = _mm_unpackhi_epi8(a0, b0);
        c1 = _mm_unpacklo_epi8(c0, d0);
        d1 = _mm_unpackhi_epi8(c0, d0);

        a0 = _mm_unpacklo_epi8(a1, b1);
        b0 = _mm_unpackhi_epi8(a1, b1);
        c0 = _mm_unpacklo_epi8(c1, d1);
        d0 = _mm_unpackhi_epi8(c1, d1);

        a1 = _mm_unpacklo_epi8(a0, b0);
        b1 = _mm_unpackhi_epi8(a0, b0);
        c1 = _mm_unpacklo_epi8(c0, d0);
        d1 = _mm_unpackhi_epi8(c0, d0);

        a0 = _mm_unpacklo_epi64(a1, c1);
        b0 = _mm_unpackhi_epi64(a1, c1);
        c0 = _mm_unpacklo_epi64(b1, d1);
        d0 = _mm_unpackhi_epi64(b1, d1);

        _mm_storeu_si128((__m128i *) &out_b[0*size + ii], a0);
        _mm_storeu_si128((__m128i *) &out_b[1*size + ii], b0);
        _mm_storeu_si128((__m128i *) &out_b[2*size + ii], c0);
        _mm_storeu_si128((__m128i *) &out_b[3*size + ii], d0);
    }
    return bshuf_trans_byte_elem_remainder(in, out, size, 4,
            size - size % 16);
}


/* Transpose bytes within elements for 64 bit elements. */
__attribute__ ((target ("sse2")))
int64_t bshuf_trans_byte_elem_SSE_64(const void* in, void* out, const size_t size) {
    
    size_t ii;
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;
    __m128i a0, b0, c0, d0, e0, f0, g0, h0;
    __m128i a1, b1, c1, d1, e1, f1, g1, h1;

    for (ii=0; ii + 15 < size; ii += 16) {
        a0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 0*16]);
        b0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 1*16]);
        c0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 2*16]);
        d0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 3*16]);
        e0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 4*16]);
        f0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 5*16]);
        g0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 6*16]);
        h0 = _mm_loadu_si128((__m128i *) &in_b[8*ii + 7*16]);

        a1 = _mm_unpacklo_epi8(a0, b0);
        b1 = _mm_unpackhi_epi8(a0, b0);
        c1 = _mm_unpacklo_epi8(c0, d0);
        d1 = _mm_unpackhi_epi8(c0, d0);
        e1 = _mm_unpacklo_epi8(e0, f0);
        f1 = _mm_unpackhi_epi8(e0, f0);
        g1 = _mm_unpacklo_epi8(g0, h0);
        h1 = _mm_unpackhi_epi8(g0, h0);

        a0 = _mm_unpacklo_epi8(a1, b1);
        b0 = _mm_unpackhi_epi8(a1, b1);
        c0 = _mm_unpacklo_epi8(c1, d1);
        d0 = _mm_unpackhi_epi8(c1, d1);
        e0 = _mm_unpacklo_epi8(e1, f1);
        f0 = _mm_unpackhi_epi8(e1, f1);
        g0 = _mm_unpacklo_epi8(g1, h1);
        h0 = _mm_unpackhi_epi8(g1, h1);

        a1 = _mm_unpacklo_epi32(a0, c0);
        b1 = _mm_unpackhi_epi32(a0, c0);
        c1 = _mm_unpacklo_epi32(b0, d0);
        d1 = _mm_unpackhi_epi32(b0, d0);
        e1 = _mm_unpacklo_epi32(e0, g0);
        f1 = _mm_unpackhi_epi32(e0, g0);
        g1 = _mm_unpacklo_epi32(f0, h0);
        h1 = _mm_unpackhi_epi32(f0, h0);

        a0 = _mm_unpacklo_epi64(a1, e1);
        b0 = _mm_unpackhi_epi64(a1, e1);
        c0 = _mm_unpacklo_epi64(b1, f1);
        d0 = _mm_unpackhi_epi64(b1, f1);
        e0 = _mm_unpacklo_epi64(c1, g1);
        f0 = _mm_unpackhi_epi64(c1, g1);
        g0 = _mm_unpacklo_epi64(d1, h1);
        h0 = _mm_unpackhi_epi64(d1, h1);

        _mm_storeu_si128((__m128i *) &out_b[0*size + ii], a0);
        _mm_storeu_si128((__m128i *) &out_b[1*size + ii], b0);
        _mm_storeu_si128((__m128i *) &out_b[2*size + ii], c0);
        _mm_storeu_si128((__m128i *) &out_b[3*size + ii], d0);
        _mm_storeu_si128((__m128i *) &out_b[4*size + ii], e0);
        _mm_storeu_si128((__m128i *) &out_b[5*size + ii], f0);
        _mm_storeu_si128((__m128i *) &out_b[6*size + ii], g0);
        _mm_storeu_si128((__m128i *) &out_b[7*size + ii], h0);
    }
    return bshuf_trans_byte_elem_remainder(in, out, size, 8,
            size - size % 16);
}


/* Transpose bytes within elements using best SSE algorithm available. */
__attribute__ ((target ("sse2")))
int64_t bshuf_trans_byte_elem_SSE(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    int64_t count;

    // Trivial cases: power of 2 bytes.
    switch (elem_size) {
        case 1:
            count = bshuf_copy(in, out, size, elem_size);
            return count;
        case 2:
            count = bshuf_trans_byte_elem_SSE_16(in, out, size);
            return count;
        case 4:
            count = bshuf_trans_byte_elem_SSE_32(in, out, size);
            return count;
        case 8:
            count = bshuf_trans_byte_elem_SSE_64(in, out, size);
            return count;
    }

    // Worst case: odd number of bytes. Turns out that this is faster for
    // (odd * 2) byte elements as well (hence % 4).
    if (elem_size % 4) {
        count = bshuf_trans_byte_elem_scal(in, out, size, elem_size);
        return count;
    }

    // Multiple of power of 2: transpose hierarchically.
    {
        size_t nchunk_elem;
        void* tmp_buf = malloc(size * elem_size);
        if (tmp_buf == NULL) return -1;

        if ((elem_size % 8) == 0) {
            nchunk_elem = elem_size / 8;
            TRANS_ELEM_TYPE(in, out, size, nchunk_elem, int64_t);
            count = bshuf_trans_byte_elem_SSE_64(out, tmp_buf,
                    size * nchunk_elem);
            bshuf_trans_elem(tmp_buf, out, 8, nchunk_elem, size);
        } else if ((elem_size % 4) == 0) {
            nchunk_elem = elem_size / 4;
            TRANS_ELEM_TYPE(in, out, size, nchunk_elem, int32_t);
            count = bshuf_trans_byte_elem_SSE_32(out, tmp_buf,
                    size * nchunk_elem);
            bshuf_trans_elem(tmp_buf, out, 4, nchunk_elem, size);
        } else {
            // Not used since scalar algorithm is faster.
            nchunk_elem = elem_size / 2;
            TRANS_ELEM_TYPE(in, out, size, nchunk_elem, int16_t);
            count = bshuf_trans_byte_elem_SSE_16(out, tmp_buf,
                    size * nchunk_elem);
            bshuf_trans_elem(tmp_buf, out, 2, nchunk_elem, size);
        }

        free(tmp_buf);
        return count;
    }
}


/* Transpose bits within bytes. */
__attribute__ ((target ("sse2")))
int64_t bshuf_trans_bit_byte_SSE(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    size_t ii, kk;
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;
    uint16_t* out_ui16;

    int64_t count;

    size_t nbyte = elem_size * size;

    CHECK_MULT_EIGHT(nbyte);

    __m128i xmm;
    int32_t bt;

    for (ii = 0; ii + 15 < nbyte; ii += 16) {
        xmm = _mm_loadu_si128((__m128i *) &in_b[ii]);
        for (kk = 0; kk < 8; kk++) {
            bt = _mm_movemask_epi8(xmm);
            xmm = _mm_slli_epi16(xmm, 1);
            out_ui16 = (uint16_t*) &out_b[((7 - kk) * nbyte + ii) / 8];
            *out_ui16 = bt;
        }
    }
    count = bshuf_trans_bit_byte_remainder(in, out, size, elem_size,
            nbyte - nbyte % 16);
    return count;
}


/* Transpose bits within elements. */
__attribute__ ((target ("sse2")))
int64_t bshuf_trans_bit_elem_SSE(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    int64_t count;

    CHECK_MULT_EIGHT(size);

    void* tmp_buf = malloc(size * elem_size);
    if (tmp_buf == NULL) return -1;

    count = bshuf_trans_byte_elem_SSE(in, out, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count = bshuf_trans_bit_byte_SSE(out, tmp_buf, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count = bshuf_trans_bitrow_eight(tmp_buf, out, size, elem_size);

    free(tmp_buf);

    return count;
}


/* For data organized into a row for each bit (8 * elem_size rows), transpose
 * the bytes. */
__attribute__ ((target ("sse2")))
int64_t bshuf_trans_byte_bitrow_SSE(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    size_t ii, jj;
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;

    CHECK_MULT_EIGHT(size);

    size_t nrows = 8 * elem_size;
    size_t nbyte_row = size / 8;

    __m128i a0, b0, c0, d0, e0, f0, g0, h0;
    __m128i a1, b1, c1, d1, e1, f1, g1, h1;
    __m128 *as, *bs, *cs, *ds, *es, *fs, *gs, *hs;

    for (ii = 0; ii + 7 < nrows; ii += 8) {
        for (jj = 0; jj + 15 < nbyte_row; jj += 16) {
            a0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 0)*nbyte_row + jj]);
            b0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 1)*nbyte_row + jj]);
            c0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 2)*nbyte_row + jj]);
            d0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 3)*nbyte_row + jj]);
            e0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 4)*nbyte_row + jj]);
            f0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 5)*nbyte_row + jj]);
            g0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 6)*nbyte_row + jj]);
            h0 = _mm_loadu_si128((__m128i *) &in_b[(ii + 7)*nbyte_row + jj]);


            a1 = _mm_unpacklo_epi8(a0, b0);
            b1 = _mm_unpacklo_epi8(c0, d0);
            c1 = _mm_unpacklo_epi8(e0, f0);
            d1 = _mm_unpacklo_epi8(g0, h0);
            e1 = _mm_unpackhi_epi8(a0, b0);
            f1 = _mm_unpackhi_epi8(c0, d0);
            g1 = _mm_unpackhi_epi8(e0, f0);
            h1 = _mm_unpackhi_epi8(g0, h0);


            a0 = _mm_unpacklo_epi16(a1, b1);
            b0 = _mm_unpacklo_epi16(c1, d1);
            c0 = _mm_unpackhi_epi16(a1, b1);
            d0 = _mm_unpackhi_epi16(c1, d1);

            e0 = _mm_unpacklo_epi16(e1, f1);
            f0 = _mm_unpacklo_epi16(g1, h1);
            g0 = _mm_unpackhi_epi16(e1, f1);
            h0 = _mm_unpackhi_epi16(g1, h1);


            a1 = _mm_unpacklo_epi32(a0, b0);
            b1 = _mm_unpackhi_epi32(a0, b0);

            c1 = _mm_unpacklo_epi32(c0, d0);
            d1 = _mm_unpackhi_epi32(c0, d0);

            e1 = _mm_unpacklo_epi32(e0, f0);
            f1 = _mm_unpackhi_epi32(e0, f0);

            g1 = _mm_unpacklo_epi32(g0, h0);
            h1 = _mm_unpackhi_epi32(g0, h0);

            // We don't have a storeh instruction for integers, so interpret
            // as a float. Have a storel (_mm_storel_epi64).
            as = (__m128 *) &a1;
            bs = (__m128 *) &b1;
            cs = (__m128 *) &c1;
            ds = (__m128 *) &d1;
            es = (__m128 *) &e1;
            fs = (__m128 *) &f1;
            gs = (__m128 *) &g1;
            hs = (__m128 *) &h1;

            _mm_storel_pi((__m64 *) &out_b[(jj + 0) * nrows + ii], *as);
            _mm_storel_pi((__m64 *) &out_b[(jj + 2) * nrows + ii], *bs);
            _mm_storel_pi((__m64 *) &out_b[(jj + 4) * nrows + ii], *cs);
            _mm_storel_pi((__m64 *) &out_b[(jj + 6) * nrows + ii], *ds);
            _mm_storel_pi((__m64 *) &out_b[(jj + 8) * nrows + ii], *es);
            _mm_storel_pi((__m64 *) &out_b[(jj + 10) * nrows + ii], *fs);
            _mm_storel_pi((__m64 *) &out_b[(jj + 12) * nrows + ii], *gs);
            _mm_storel_pi((__m64 *) &out_b[(jj + 14) * nrows + ii], *hs);

            _mm_storeh_pi((__m64 *) &out_b[(jj + 1) * nrows + ii], *as);
            _mm_storeh_pi((__m64 *) &out_b[(jj + 3) * nrows + ii], *bs);
            _mm_storeh_pi((__m64 *) &out_b[(jj + 5) * nrows + ii], *cs);
            _mm_storeh_pi((__m64 *) &out_b[(jj + 7) * nrows + ii], *ds);
            _mm_storeh_pi((__m64 *) &out_b[(jj + 9) * nrows + ii], *es);
            _mm_storeh_pi((__m64 *) &out_b[(jj + 11) * nrows + ii], *fs);
            _mm_storeh_pi((__m64 *) &out_b[(jj + 13) * nrows + ii], *gs);
            _mm_storeh_pi((__m64 *) &out_b[(jj + 15) * nrows + ii], *hs);
        }
        for (jj = nbyte_row - nbyte_row % 16; jj < nbyte_row; jj ++) {
            out_b[jj * nrows + ii + 0] = in_b[(ii + 0)*nbyte_row + jj];
            out_b[jj * nrows + ii + 1] = in_b[(ii + 1)*nbyte_row + jj];
            out_b[jj * nrows + ii + 2] = in_b[(ii + 2)*nbyte_row + jj];
            out_b[jj * nrows + ii + 3] = in_b[(ii + 3)*nbyte_row + jj];
            out_b[jj * nrows + ii + 4] = in_b[(ii + 4)*nbyte_row + jj];
            out_b[jj * nrows + ii + 5] = in_b[(ii + 5)*nbyte_row + jj];
            out_b[jj * nrows + ii + 6] = in_b[(ii + 6)*nbyte_row + jj];
            out_b[jj * nrows + ii + 7] = in_b[(ii + 7)*nbyte_row + jj];
        }
    }
    return size * elem_size;
}


/* Shuffle bits within the bytes of eight element blocks. */
__attribute__ ((target ("sse2")))
int64_t bshuf_shuffle_bit_eightelem_SSE(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    CHECK_MULT_EIGHT(size);

    // With a bit of care, this could be written such that such that it is
    // in_buf = out_buf safe.
    const char* in_b = (const char*) in;
    uint16_t* out_ui16 = (uint16_t*) out;

    size_t ii, jj, kk;
    size_t nbyte = elem_size * size;

    __m128i xmm;
    int32_t bt;

    if (elem_size % 2) {
        bshuf_shuffle_bit_eightelem_scal(in, out, size, elem_size);
    } else {
        for (ii = 0; ii + 8 * elem_size - 1 < nbyte;
                ii += 8 * elem_size) {
            for (jj = 0; jj + 15 < 8 * elem_size; jj += 16) {
                xmm = _mm_loadu_si128((__m128i *) &in_b[ii + jj]);
                for (kk = 0; kk < 8; kk++) {
                    bt = _mm_movemask_epi8(xmm);
                    xmm = _mm_slli_epi16(xmm, 1);
                    size_t ind = (ii + jj / 8 + (7 - kk) * elem_size);
                    out_ui16[ind / 2] = bt;
                }
            }
        }
    }
    return size * elem_size;
}


/* Untranspose bits within elements. */
__attribute__ ((target ("sse2")))
int64_t bshuf_untrans_bit_elem_SSE(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    int64_t count;

    CHECK_MULT_EIGHT(size);

    void* tmp_buf = malloc(size * elem_size);
    if (tmp_buf == NULL) return -1;

    count = bshuf_trans_byte_bitrow_SSE(in, tmp_buf, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count =  bshuf_shuffle_bit_eightelem_SSE(tmp_buf, out, size, elem_size);

    free(tmp_buf);

    return count;
}


/* ---- Code that requires AVX2. Intel Haswell (2013) and later. ---- */

/* ---- Worker code that uses AVX2 ----
 *
 * The following code makes use of the AVX2 instruction set and specialized
 * 32 byte registers. The AVX2 instructions are present on newer x86
 * processors. The first Intel processor microarchitecture supporting AVX2 was
 * Haswell (2013).
 *
 */

/* Transpose bits within bytes. */
__attribute__ ((target ("avx2")))
int64_t bshuf_trans_bit_byte_AVX(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    size_t ii, kk;
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;
    int32_t* out_i32;

    size_t nbyte = elem_size * size;

    int64_t count;

    __m256i ymm;
    int32_t bt;

    for (ii = 0; ii + 31 < nbyte; ii += 32) {
        ymm = _mm256_loadu_si256((__m256i *) &in_b[ii]);
        for (kk = 0; kk < 8; kk++) {
            bt = _mm256_movemask_epi8(ymm);
            ymm = _mm256_slli_epi16(ymm, 1);
            out_i32 = (int32_t*) &out_b[((7 - kk) * nbyte + ii) / 8];
            *out_i32 = bt;
        }
    }
    count = bshuf_trans_bit_byte_remainder(in, out, size, elem_size,
            nbyte - nbyte % 32);
    return count;
}


/* Transpose bits within elements. */
__attribute__ ((target ("avx2")))
int64_t bshuf_trans_bit_elem_AVX(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    int64_t count;

    CHECK_MULT_EIGHT(size);

    void* tmp_buf = malloc(size * elem_size);
    if (tmp_buf == NULL) return -1;

    count = bshuf_trans_byte_elem_SSE(in, out, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count = bshuf_trans_bit_byte_AVX(out, tmp_buf, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count = bshuf_trans_bitrow_eight(tmp_buf, out, size, elem_size);

    free(tmp_buf);

    return count;
}


/* For data organized into a row for each bit (8 * elem_size rows), transpose
 * the bytes. */
__attribute__ ((target ("avx2")))
int64_t bshuf_trans_byte_bitrow_AVX(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    size_t hh, ii, jj, kk, mm;
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;

    CHECK_MULT_EIGHT(size);

    size_t nrows = 8 * elem_size;
    size_t nbyte_row = size / 8;

    if (elem_size % 4) return bshuf_trans_byte_bitrow_SSE(in, out, size,
            elem_size);

    __m256i ymm_0[8];
    __m256i ymm_1[8];
    __m256i ymm_storeage[8][4];

    for (jj = 0; jj + 31 < nbyte_row; jj += 32) {
        for (ii = 0; ii + 3 < elem_size; ii += 4) {
            for (hh = 0; hh < 4; hh ++) {

                for (kk = 0; kk < 8; kk ++){
                    ymm_0[kk] = _mm256_loadu_si256((__m256i *) &in_b[
                            (ii * 8 + hh * 8 + kk) * nbyte_row + jj]);
                }

                for (kk = 0; kk < 4; kk ++){
                    ymm_1[kk] = _mm256_unpacklo_epi8(ymm_0[kk * 2],
                            ymm_0[kk * 2 + 1]);
                    ymm_1[kk + 4] = _mm256_unpackhi_epi8(ymm_0[kk * 2],
                            ymm_0[kk * 2 + 1]);
                }

                for (kk = 0; kk < 2; kk ++){
                    for (mm = 0; mm < 2; mm ++){
                        ymm_0[kk * 4 + mm] = _mm256_unpacklo_epi16(
                                ymm_1[kk * 4 + mm * 2],
                                ymm_1[kk * 4 + mm * 2 + 1]);
                        ymm_0[kk * 4 + mm + 2] = _mm256_unpackhi_epi16(
                                ymm_1[kk * 4 + mm * 2],
                                ymm_1[kk * 4 + mm * 2 + 1]);
                    }
                }

                for (kk = 0; kk < 4; kk ++){
                    ymm_1[kk * 2] = _mm256_unpacklo_epi32(ymm_0[kk * 2],
                            ymm_0[kk * 2 + 1]);
                    ymm_1[kk * 2 + 1] = _mm256_unpackhi_epi32(ymm_0[kk * 2],
                            ymm_0[kk * 2 + 1]);
                }

                for (kk = 0; kk < 8; kk ++){
                    ymm_storeage[kk][hh] = ymm_1[kk];
                }
            }

            for (mm = 0; mm < 8; mm ++) {

                for (kk = 0; kk < 4; kk ++){
                    ymm_0[kk] = ymm_storeage[mm][kk];
                }

                ymm_1[0] = _mm256_unpacklo_epi64(ymm_0[0], ymm_0[1]);
                ymm_1[1] = _mm256_unpacklo_epi64(ymm_0[2], ymm_0[3]);
                ymm_1[2] = _mm256_unpackhi_epi64(ymm_0[0], ymm_0[1]);
                ymm_1[3] = _mm256_unpackhi_epi64(ymm_0[2], ymm_0[3]);

                ymm_0[0] = _mm256_permute2x128_si256(ymm_1[0], ymm_1[1], 32);
                ymm_0[1] = _mm256_permute2x128_si256(ymm_1[2], ymm_1[3], 32);
                ymm_0[2] = _mm256_permute2x128_si256(ymm_1[0], ymm_1[1], 49);
                ymm_0[3] = _mm256_permute2x128_si256(ymm_1[2], ymm_1[3], 49);

                _mm256_storeu_si256((__m256i *) &out_b[
                        (jj + mm * 2 + 0 * 16) * nrows + ii * 8], ymm_0[0]);
                _mm256_storeu_si256((__m256i *) &out_b[
                        (jj + mm * 2 + 0 * 16 + 1) * nrows + ii * 8], ymm_0[1]);
                _mm256_storeu_si256((__m256i *) &out_b[
                        (jj + mm * 2 + 1 * 16) * nrows + ii * 8], ymm_0[2]);
                _mm256_storeu_si256((__m256i *) &out_b[
                        (jj + mm * 2 + 1 * 16 + 1) * nrows + ii * 8], ymm_0[3]);
            }
        }
    }
    for (ii = 0; ii < nrows; ii ++ ) {
        for (jj = nbyte_row - nbyte_row % 32; jj < nbyte_row; jj ++) {
            out_b[jj * nrows + ii] = in_b[ii * nbyte_row + jj];
        }
    }
    return size * elem_size;
}


/* Shuffle bits within the bytes of eight element blocks. */
__attribute__ ((target ("avx2")))
int64_t bshuf_shuffle_bit_eightelem_AVX(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    CHECK_MULT_EIGHT(size);

    // With a bit of care, this could be written such that such that it is
    // in_buf = out_buf safe.
    const char* in_b = (const char*) in;
    char* out_b = (char*) out;

    size_t ii, jj, kk;
    size_t nbyte = elem_size * size;

    __m256i ymm;
    int32_t bt;

    if (elem_size % 4) {
        return bshuf_shuffle_bit_eightelem_SSE(in, out, size, elem_size);
    } else {
        for (jj = 0; jj + 31 < 8 * elem_size; jj += 32) {
            for (ii = 0; ii + 8 * elem_size - 1 < nbyte;
                    ii += 8 * elem_size) {
                ymm = _mm256_loadu_si256((__m256i *) &in_b[ii + jj]);
                for (kk = 0; kk < 8; kk++) {
                    bt = _mm256_movemask_epi8(ymm);
                    ymm = _mm256_slli_epi16(ymm, 1);
                    size_t ind = (ii + jj / 8 + (7 - kk) * elem_size);
                    * (int32_t *) &out_b[ind] = bt;
                }
            }
        }
    }
    return size * elem_size;
}


/* Untranspose bits within elements. */
__attribute__ ((target ("avx2")))
int64_t bshuf_untrans_bit_elem_AVX(const void* in, void* out, const size_t size,
         const size_t elem_size) {
    
    int64_t count;

    CHECK_MULT_EIGHT(size);

    void* tmp_buf = malloc(size * elem_size);
    if (tmp_buf == NULL) return -1;

    count = bshuf_trans_byte_bitrow_AVX(in, tmp_buf, size, elem_size);
    CHECK_ERR_FREE(count, tmp_buf);
    count =  bshuf_shuffle_bit_eightelem_AVX(tmp_buf, out, size, elem_size);

    free(tmp_buf);
    return count;
}


/* ---- Drivers selecting best instruction set at run time. ---- */

int64_t bshuf_trans_bit_elem(const void* in, void* out, const size_t size, const size_t elem_size) {
    unsigned int eax, ebx, ecx, edx;
	signed char has_sse2 = 0;

    /* Collect CPU features */
    if (!__get_cpuid (1, &eax, &ebx, &ecx, &edx)){
        return bshuf_trans_bit_elem_scal(in, out, size, elem_size);
    }

	has_sse2 = ((edx & bit_SSE2) != 0);
    if (__get_cpuid_max (0, NULL) < 7){
        return bshuf_trans_bit_elem_scal(in, out, size, elem_size);
    }

	__cpuid_count (7, 0, eax, ebx, ecx, edx);
    
	signed char has_avx2 = 0;
	has_avx2 = ((ebx & bit_AVX2) != 0);
	if (has_avx2){
	    return bshuf_trans_bit_elem_AVX(in, out, size, elem_size);
    }

    /* Pick SSE2 version */
    if (has_sse2){
        return bshuf_trans_bit_elem_SSE(in, out, size, elem_size);
    }
    /* Fallback to default implementation */
    return bshuf_trans_bit_elem_scal(in, out, size, elem_size);
}

int64_t bshuf_untrans_bit_elem(const void* in, void* out, const size_t size, const size_t elem_size) {
    unsigned int eax, ebx, ecx, edx;
	signed char has_sse2 = 0;

    /* Collect CPU features */
    if (!__get_cpuid (1, &eax, &ebx, &ecx, &edx)){
        return bshuf_untrans_bit_elem_scal(in, out, size, elem_size);
    }

	has_sse2 = ((edx & bit_SSE2) != 0);
    if (__get_cpuid_max (0, NULL) < 7){
        return bshuf_untrans_bit_elem_scal(in, out, size, elem_size);
    }

	__cpuid_count (7, 0, eax, ebx, ecx, edx);
    
	signed char has_avx2 = 0;
	has_avx2 = ((ebx & bit_AVX2) != 0);
	if (has_avx2){
	    return bshuf_untrans_bit_elem_AVX(in, out, size, elem_size);
    }

    /* Pick SSE2 version */
    if (has_sse2){
        return bshuf_untrans_bit_elem_SSE(in, out, size, elem_size);
    }
    /* Fallback to default implementation */
    return bshuf_untrans_bit_elem_scal(in, out, size, elem_size);
}

#else // x86_64

int64_t bshuf_trans_bit_elem(const void* in, void* out, const size_t size, const size_t elem_size){
    return bshuf_trans_bit_elem_scal(in, out, size, elem_size);
}

int64_t bshuf_untrans_bit_elem(const void* in, void* out, const size_t size, const size_t elem_size){
    return bshuf_untrans_bit_elem_scal(in, out, size, elem_size);
}

#endif


/* ---- Wrappers for implementing blocking ---- */

/* Wrap a function for processing a single block to process an entire buffer in
 * parallel. */
int64_t bshuf_blocked_wrap_fun(bshufBlockFunDef fun, const void* in, void* out, \
        const size_t size, const size_t elem_size, size_t block_size) {

    omp_size_t ii = 0;
    int64_t err = 0;
    int64_t count, cum_count=0;
    size_t last_block_size;
    size_t leftover_bytes;
    size_t this_iter;
    char *last_in;
    char *last_out;


    ioc_chain C;
    ioc_init(&C, in, out);


    if (block_size == 0) {
        block_size = bshuf_default_block_size(elem_size);
    }
    if (block_size % BSHUF_BLOCKED_MULT) return -81;

    for (ii = 0; ii < (omp_size_t)( size / block_size ); ii ++) {
        count = fun(&C, block_size, elem_size);
        if (count < 0) err = count;
        cum_count += count;
    }

    last_block_size = size % block_size;
    last_block_size = last_block_size - last_block_size % BSHUF_BLOCKED_MULT;
    if (last_block_size) {
        count = fun(&C, last_block_size, elem_size);
        if (count < 0) err = count;
        cum_count += count;
    }

    if (err < 0) return err;

    leftover_bytes = size % BSHUF_BLOCKED_MULT * elem_size;
    last_in = (char *) ioc_get_in(&C, &this_iter);
    ioc_set_next_in(&C, &this_iter, (void *) (last_in + leftover_bytes));
    last_out = (char *) ioc_get_out(&C, &this_iter);
    ioc_set_next_out(&C, &this_iter, (void *) (last_out + leftover_bytes));

    memcpy(last_out, last_in, leftover_bytes);

    ioc_destroy(&C);

    return cum_count + leftover_bytes;
}


/* Write a 32 bit unsigned integer to a buffer in big endian order. */
void bshuf_write_uint32_BE(void* buf, uint32_t num) {
    int ii;
    uint8_t* b = (uint8_t*) buf;
    uint32_t pow28 = 1 << 8;
    for (ii = 3; ii >= 0; ii--) {
        b[ii] = num % pow28;
        num = num / pow28;
    }
}


/* Read a 32 bit unsigned integer from a buffer big endian order. */
uint32_t bshuf_read_uint32_BE(const void* buf) {
    int ii;
    uint8_t* b = (uint8_t*) buf;
    uint32_t num = 0, pow28 = 1 << 8, cp = 1;
    for (ii = 3; ii >= 0; ii--) {
        num += b[ii] * cp;
        cp *= pow28;
    }
    return num;
}


/* ---- Public functions ----
 *
 * See header file for description and usage.
 *
 */

size_t bshuf_default_block_size(const size_t elem_size) {
    // This function needs to be absolutely stable between versions.
    // Otherwise encoded data will not be decodable.

    size_t block_size = BSHUF_TARGET_BLOCK_SIZE_B / elem_size;
    // Ensure it is a required multiple.
    block_size = (block_size / BSHUF_BLOCKED_MULT) * BSHUF_BLOCKED_MULT;
    return MAX(block_size, BSHUF_MIN_RECOMMEND_BLOCK);
}


#undef TRANS_BIT_8X8
#undef TRANS_ELEM_TYPE
#undef MAX
#undef CHECK_MULT_EIGHT
#undef CHECK_ERR_FREE
