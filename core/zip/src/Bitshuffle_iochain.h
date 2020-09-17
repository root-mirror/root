/*
 * IOchain - Distribute a chain of dependant IO events amoung threads.
 *
 * This file is part of Bitshuffle
 * Author: Kiyoshi Masui <kiyo@physics.ubc.ca>
 * Website: http://www.github.com/kiyo-masui/bitshuffle
 * Created: 2014
 *
 * See LICENSE file for details about copyright and rights to use.
 *
 *
 * Header File
 *
 * Similar in concept to a queue. Each task includes reading an input
 * and writing output, but the location of the input/output (the pointers)
 * depend on the previous item in the chain.
 *
 * This is designed for parallelizing blocked compression/decompression IO,
 * where the destination of a compressed block depends on the compressed size
 * of all previous blocks.
 *
 * Implemented with OpenMP locks.
 *
 *
 * Usage
 * -----
 *  - Call `ioc_init` in serial block.
 *  - Each thread should create a local variable *size_t this_iter* and
 *    pass its address to all function calls. Its value will be set
 *    inside the functions and is used to identify the thread.
 *  - Each thread must call each of the `ioc_get*` and `ioc_set*` methods
 *    exactly once per iteration, starting with `ioc_get_in` and ending
 *    with `ioc_set_next_out`.
 *  - The order (`ioc_get_in`, `ioc_set_next_in`, *work*, `ioc_get_out`,
 *    `ioc_set_next_out`, *work*) is most efficient.
 *  - Have each thread call `ioc_end_pop`.
 *  - `ioc_get_in` is blocked until the previous entry's
 *    `ioc_set_next_in` is called.
 *  - `ioc_get_out` is blocked until the previous entry's
 *    `ioc_set_next_out` is called.
 *  - There are no blocks on the very first iteration.
 *  - Call `ioc_destroy` in serial block.
 *  - Safe for num_threads >= IOC_SIZE (but less efficient).
 *
 */


#ifndef IOCHAIN_H
#define IOCHAIN_H


#include <stdlib.h>


#define IOC_SIZE 33


typedef struct ioc_ptr_and_lock {
    void *ptr;
} ptr_and_lock;

typedef struct ioc_const_ptr_and_lock {
    const void *ptr;
} const_ptr_and_lock;


typedef struct ioc_chain {
    size_t next;
    const_ptr_and_lock in_pl[IOC_SIZE];
    ptr_and_lock out_pl[IOC_SIZE];
} ioc_chain;


void ioc_init(ioc_chain *C, const void *in_ptr_0, void *out_ptr_0);
void ioc_destroy(ioc_chain *C);
const void * ioc_get_in(ioc_chain *C, size_t *this_iter);
void ioc_set_next_in(ioc_chain *C, size_t* this_iter, void* in_ptr);
void * ioc_get_out(ioc_chain *C, size_t *this_iter);
void ioc_set_next_out(ioc_chain *C, size_t *this_iter, void* out_ptr);

#endif  // IOCHAIN_H
