/*
 * bn256-20080525/fp12.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef FP12_H
#define FP12_H

#include "fp6e.h"

#ifdef BENCH
unsigned long long multp12cycles; unsigned long long nummultp12;
unsigned long long sqp12cycles; unsigned long long numsqp12;
unsigned long long invp12cycles; unsigned long long numinvp12;
#endif

extern fp6e_t tau; // F_{p^{12}} is constructed as F_{p^6}[Z]/(Z^2 - tau) F_{p^6}[Z]
extern fp2e_t zpminus1; // Z^{p-1}, lies in F_{p^2}
extern fp2e_t zpminus1inv; // Z^{p-1}, lies in F_{p^2}

void fp12_init(void);

void fp12_clear(void);

#endif // ifdef FP12_H
