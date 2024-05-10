/*
 * bn256-20080525/curve.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef CURVE_H
#define CURVE_H

#include <gmp.h>

extern mpz_t b; /* parameter b in the curve equation y^2 = x^3 + b */
extern mpz_t n; /* order of the curve */
extern mpz_t x;  
extern mpz_t trace; /* trace of Frobenius */
extern mpz_t chi; /* p^12 / n */
extern mpz_t chicomp; 
extern mpz_t looplength_eta; 

/* Initialize Curve paramesters */
void curve_init(void);

/* Free memory used by curve parameters */
void curve_clear(void);

#endif /* ifdef CURVE_H */
