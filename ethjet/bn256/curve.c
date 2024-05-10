/*
 * bn256-20080525/curve.c 
 * Peter Schwabe
 * Public domain
*/

#include "curve.h"
#include "parameters.h"

mpz_t b; /* parameter b in the curve equation y^2 = x^3 + b */
mpz_t n; /* order of the curve */
mpz_t x;  
mpz_t trace; /* trace of Frobenius */
mpz_t chi; /* p^12 / n */
mpz_t chicomp; 
mpz_t looplength_eta;

void curve_init()
{
	/* Curve parameters */
	mpz_init_set_str(x, BN_X, 10);
	mpz_init_set_str(n, BN_N, 10);
	mpz_init_set_str(trace, BN_TRACE,10);
	mpz_init_set_str(chi, BN_CHI, 10); // (p^k - 1) / n
	mpz_init_set_str(chicomp, BN_CHICOMP, 10);
	mpz_init_set_str(looplength_eta, BN_LOOPLENGTH_ETA, 10);
	mpz_init_set_str(b, BN_B, 10);
}

void curve_clear()
{
	mpz_clear(n);
	mpz_clear(trace);
	mpz_clear(chicomp);
	mpz_clear(chi);
	mpz_clear(x);
	mpz_clear(b);
	mpz_clear(looplength_eta);
}
