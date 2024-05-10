/*
 * bn256-20080525/fp.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef FP_H
#define FP_H

#include <gmp.h>

#include "fpe.h"

#ifdef BENCH
unsigned long long int multpcycles; unsigned long long int nummultp;
unsigned long long int sqpcycles; unsigned long long int numsqp;
unsigned long long invpcycles; unsigned long long numinvp;
#endif


extern mpz_t p;
extern unsigned long long p_inv; // -p^{-1} mod 2^{GMP_LIMB_BITS} used in Montgomery reduction

void fp_init(void);

void fp_clear(void);

#endif // ifdef FP_H
