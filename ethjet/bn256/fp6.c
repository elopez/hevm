/*
 * bn256-20080525/fp6.c 
 * Peter Schwabe
 * Public domain
*/

#include "parameters.h"
#include "fp6.h"

fp2e_t xi; // constant coefficient in the irreducible polynomial Y^3 - xi, used to construct F_{p^6} as extension of F_{p^12}
fp2e_t ypminus1; // Y^{p-1} lies in F_{p^2}

fpe_t zeta; // Third root of unity in F_p fulfilling Z^{p^2} = -zeta * Z

// Multiples of powers of xi needed for cometa-pairing computation
fp2e_t xi2; // xi^2
fp2e_t _1o27xi3; // 1/27 xi^3
fp2e_t _1o3xi3; // 1/3 xi^3
fp2e_t _1o3xi; // 1/3 xi
fpe_t _1o3modp; // 1/3 \in \F_p

// Two constants needed for the cometa-pairing computation
fpe_t cometa_c0_const;
fpe_t cometa_c1_const;

void fp6_init()
{
	fp2e_set_str(xi, BN_XI);
	fp2e_set_str(ypminus1, BN_YPMINUS1);
	fpe_set_str(zeta, BN_ZETA);
	fp2e_set_str(xi2, BN_XI2);
	fp2e_set_str(_1o27xi3, BN_1O27XI3);
	fp2e_set_str(_1o3xi3, BN_1O3XI3);
	fp2e_set_str(_1o3xi, BN_1O3XI);
	fpe_set_str(_1o3modp, BN_1O3MODP); // 1/3 \in \F_p
	
	fpe_set_str(cometa_c0_const, BN_COMETA_C0_CONST);
	fpe_set_str(cometa_c1_const, BN_COMETA_C1_CONST);
}

void fp6_clear()
{
	;
}
