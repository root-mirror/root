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
 */

#include <stdlib.h>
#include "Bitshuffle_iochain.h"


void ioc_init(ioc_chain *C, const void *in_ptr_0, void *out_ptr_0) {
    C->next = 0;
    C->in_pl[0].ptr = in_ptr_0;
    C->out_pl[0].ptr = out_ptr_0;
}


void ioc_destroy(ioc_chain *C) {
    (void)C;
}


const void * ioc_get_in(ioc_chain *C, size_t *this_iter) {
    *this_iter = C->next;
    C->next ++;
    return C->in_pl[*this_iter % IOC_SIZE].ptr;
}


void ioc_set_next_in(ioc_chain *C, size_t* this_iter, void* in_ptr) {
    C->in_pl[(*this_iter + 1) % IOC_SIZE].ptr = in_ptr;
}


void * ioc_get_out(ioc_chain *C, size_t *this_iter) {
    void *out_ptr = C->out_pl[*this_iter % IOC_SIZE].ptr;
    return out_ptr;
}


void ioc_set_next_out(ioc_chain *C, size_t *this_iter, void* out_ptr) {
    C->out_pl[(*this_iter + 1) % IOC_SIZE].ptr = out_ptr;
}
