/*
 * bn256-20080525/fp2.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef FP2_H
#define FP2_H

#include <gmp.h>
#include "fp.h"
#include "parameters.h"

#ifdef BENCH
unsigned long long multp2cycles; unsigned long long nummultp2;
unsigned long long sqp2cycles; unsigned long long numsqp2;
unsigned long long invp2cycles; unsigned long long numinvp2;
#endif

void fp2_init(void);
void fp2_clear(void);

#endif // ifdef FP2_H
