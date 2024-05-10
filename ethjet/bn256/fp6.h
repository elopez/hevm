/*
 * bn256-20080525/fp6.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef FP6_H
#define FP6_H

#include "fp2e.h"

extern fp2e_t xi; // constant coefficient in the irreducible polynomial Y^3 - xi, used to construct F_{p^6} as extension of F_{p^12}
extern fp2e_t ypminus1; // Y^{p-1} lies in F_{p^2}

extern fpe_t zeta; // Third root of unity in F_p fulfilling Z^{p^2} = -zeta * Z

// Multiples of powers of xi needed for cometa-pairing computation
extern fp2e_t xi2; // xi^2
extern fp2e_t _1o27xi3; // 1/27 xi^3
extern fp2e_t _1o3xi3; // 1/3 xi^3
extern fp2e_t _1o3xi; // 1/3 xi
extern fpe_t _1o3modp; // 1/3 \in \F_p

// Two constants needed for the cometa-pairing computation
extern fpe_t cometa_c0_const;
extern fpe_t cometa_c1_const;

void fp6_init(void);
void fp6_clear(void);

#endif // ifdef FP6_H
